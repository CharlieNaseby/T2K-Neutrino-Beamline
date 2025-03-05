import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import re
import copy as cp
import sys
import ast
#sys.path.append('../')
#import create_beamline
#matplotlib.use('qtagg')
pd.set_option('display.max_rows', 500000)
#pd.set_option('display.max_columns', 500000)


#known constants

start_of_bpd1 = 17400.
end_of_bpd2 = 24598.
bpd1_angle = 0.033547
bpd2_angle = 0.033165
ssem1s = 1769.-85.
proton_momentum = 30.924 # momentum for a 30GeV KE proton 
vacuum_pressure = 1e-4 #vacuum pressure in bar

#print_tunnel=False
#print_physics=True
#sample_all=True
#sample_ssem=False
#sample_entry=False
#beam_from_file = False
#beam_halo = False
#enable_blms = False
#geometry = True
#misalignments = False
#print_vacuum=True
#bias_physics=True

#fit configuration
print_tunnel=False
print_physics=False
sample_all=False
sample_ssem=True
sample_entry=False  #####WARNING MUST BE FALSE WHEN FITTING OTHERWISE ENTRY WILL BE TREATED AS SSEM1!!!
beam_from_file = False
beam_halo = False
enable_blms = False
geometry = False
misalignments=True
bias_physics=False
print_vacuum=False
merge_drifts=False

generate_primaries=False

if(generate_primaries):
    sample_entry = True
    sample_ssem = False
    sample_all = False
    print_physics = False
    print_tunnel = False
    beam_from_file = False
    beam_halo = False
    enable_blms = False
    bias_physics=False

def get_rotation_matrix(a, b):
    #first normalise the vectors just to be safe
    a /= np.linalg.norm(a)
    b /= np.linalg.norm(b)
   
    #function to find the rotation matrix that maps vector a to b
    #https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    #NOTE THERE IS AN ADDITIONAL DEGREE OF FREEDOM - ROTATION ABOUT THE AXIS OF THE VECTOR this must be dealt with carefully
    #sadly we only have two points on each component (some we only have one!) so this is a perennial problem
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = a.dot(b)
    vx = np.array([[0., -v[2], v[1]],
                    [v[2], 0., -v[0]],
                    [-v[1], v[0], 0.]])
    R = np.identity(3) + vx + vx.dot(vx)*1./(1.+c)
    return R

def is_empty(obj):
    if isinstance(obj, np.ndarray):  # If it's a NumPy array
        return obj.size == 0  # Check if the array is empty using its size
    elif isinstance(obj, list):  # If it's a Python list
        return len(obj) == 0  # Check if the list is empty using its length
    else:
        raise ValueError("Unsupported type. Only lists and numpy arrays are supported.")
        


def rotmat2d(theta):
    return np.array([[np.cos(theta), np.sin(theta)],
                    [-np.sin(theta), np.cos(theta)]])

def normvec(a):
    retval = np.array(a)
    retval /= np.linalg.norm(retval)
    return retval

def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

def get_intersection_len(p1, v1, p2, v2):
    A = np.array([[v1[0], -v2[0]], [v1[1], -v2[1]]])
    B = p2 - p1
    a, b = np.linalg.solve(A, B)
    return a
    
    
#    umerator = p2[1]-p1[1]+v2[1]*(p1[0]-p2[0])/v2[0]
#    denominator = v1[1]+v1[0]*v2[1]/v2[0]
#    return numerator/denominator


def string_to_list(s):
    s = s.replace(' ', ',')
    s = s.replace('][', '],[')
    elements = ast.literal_eval(s)
    return elements

class survey_element:
    def __init__(self, nom_normvec):
        self.nominal_normvec = cp.deepcopy(nom_normvec)
        self.s_nom = 0.
        self.s_surv = 0.
        self.segment = []
        self.segment_id = -1.
        self.nominal_center = []
        self.survey_points = []
        self.survey_center = []
        self.survey_diff = []
        self.nominal_offset = []
        self.name = []

    def rotate_survey(self, axis, theta):
        c = np.cos(theta)
        s = np.sin(theta)
        if(axis == 'xy'):
            R = np.array([[c, -s, 0.],[s, c, 0.], [0., 0., 1.]])
        else:
            raise Exception('axis not valid')
        if(not is_empty(self.survey_points)):
            self.survey_points = R.dot(self.survey_points.T).T
        if(not is_empty(self.survey_center)):
            self.survey_center = R.dot(self.survey_center.T).T
    def shift_survey(self, change):
        self.survey_center += change
        self.survey_points += change
    
    def binary_search(self, func, R, horisontal, l, r):
        precision = 1e-8
        val = 99999
        left = [l, func(l, R, horisontal)]
        right = [r, func(r, R, horisontal)]
        niter = 0
        while(np.abs(val) > precision and niter < 1000):
            niter += 1
            testpoint = (left[0]+right[0])/2.
            val = func(testpoint, R, horisontal)
#            print(f'testing point {testpoint} returned value {val}')
            if(val<0):
                if(left[1]<0):
                    left = [testpoint, val]
                else:
                    right = [testpoint, val]
            elif(val>0):
                if(left[1]>0):
                    left = [testpoint, val]
                else:
                    right = [testpoint, val]
        if(np.abs(right[1]) < np.abs(left[1])):
            return right[0]
        else:
            return left[0]

    def calculate_survey_center(self, s):

        survey_point_diff = self.survey_points[1] - self.survey_points[0]
        nominal_surv_point_diff = -self.nominal_offset[1] + self.nominal_offset[0] #these numbers are center to point vector so sort of backwards from the survey

        R = get_rotation_matrix(nominal_surv_point_diff, survey_point_diff) #rotation to take the points from the idealised space to the survey space

        p1c = R.dot(-self.nominal_offset[0]) #direction to go from survey point 1 to the center
        #there is an extra DOF, rotation about the line connecting the two survey points traces out a cone of valid p1c vectors
        #how far along the magnet direction does the p1c line go?
        len_along_magdir = p1c.dot(normvec(-survey_point_diff))
#        print(f'len along magdir = {len_along_magdir}')
        circle_radius = np.linalg.norm(np.cross(p1c, normvec(survey_point_diff)))

        #can now define a circle around survey_point_diff at position len_along_magdir with radius circle_radius
        R2 = get_rotation_matrix(np.array([0.0, 0.0, 1.0]), survey_point_diff)
        angles = np.linspace(0, 2*np.pi, 100)
        xycircle = np.array([np.cos(angles), np.sin(angles), np.zeros(angles.size)])
        rotated_circle = R2.dot(xycircle)
        vertical = np.array([0.0, 0.0, 1.0])
        horisontal = normvec(np.cross(vertical, normvec(survey_point_diff)))
        circle = lambda angle, R2, horisontal: horisontal.dot(R2.dot(np.array([np.cos(angle), np.sin(angle), 0.0])))

        best_angle = self.binary_search(circle, R2, horisontal, -np.pi/2., np.pi/2.)
#        print(f'best angle {best_angle}')
        #get resulting center point
        self.survey_center = self.survey_points[0] + len_along_magdir*normvec(survey_point_diff) + circle_radius*R2.dot(np.array([np.cos(best_angle), np.sin(best_angle), 0.0]))
#        print(f'survey center {self.survey_center}')
#        print(f'for survey points {self.survey_points}')

        #alternate method as a cross-check
        center = self.survey_points[0]
        center -= self.nominal_offset[0][0]*normvec(survey_point_diff) #move along vector connecting survey points by amount nominal diff
        center -= self.nominal_offset[0][1]*horisontal
        center -= self.nominal_offset[0][2]*vertical
#        print(f'alternative method gets center of {center}')


class nominalBeamline:
    svectors = []
    sstart = []
    slen = []
    bline = []
    offsets = []
    mag_objs = []
    def __init__(self): #a bunch of vectors extraced from the .dwg files provided by Fujii-san
        self.nuin_to_bpd1 = np.array([[27711.320296225927,10246.557442498312,0.0], [8826.078822922398,10998.629426709593,0.0]])
        self.startpoint = self.nuin_to_bpd1[0,:]
        self.bpd1_to_bpd2 = np.array([[8826.078822922398,10998.629426709593,0.0], [4625.740813728273,11025.057932632815,0.0]])
        self.bpd2_to_arc = np.array([[4625.740813728273,11025.057932632815,0.0],[-26953.86986321407,10165.298827699111,0.0]])
        self.mag_objs = {"BPV1":survey_element(np.array([[24832.38910462714,9576.830981585355,0.0], [24885.262495664367,10904.532021424293,0.0]])),
                         "BPH2": survey_element(np.array([[20234.464068271784,10976.996762330913,0.0], [20189.0617097797,9836.900439275993,0.0]])),
                         "QPQ1": survey_element(np.array([[17308.21253825885,11632.058913716024,0.0], [17231.215549930905,9698.591442102286,0.0]])),
                         "QPQ2": survey_element(np.array([[12997.315158351666,11444.001473098717,0.0], [12947.273757007451,10187.414380095572,0.0]])),
                         "BPD1": survey_element(np.array([[8849.395131864345,12038.357324908586,0.0], [8802.84416590732,10018.893780319399,0.0]])),
                         "BPD2": survey_element(np.array([[4614.8592335086605,12065.001003860916,0.0], [4635.994610472831,10045.111577135578,0.0]])),
                         "QPQ3": survey_element(np.array([[-5990.328971827361,11463.303398540442,0.0], [-5955.518885550551,10184.702260589635,0.0]])),
                         "BPV2": survey_element(np.array([[-10385.778916663812,11237.333989977218,0.0], [-10351.093606239794,9963.315964633332,0.0]])),
                         "QPQ4": survey_element(np.array([[-13482.390984855087,11437.10340747186,0.0], [-13436.376745104855,9746.96512861545,0.0]])),
                         "BPH3": survey_element(np.array([[-21036.911835919527,11136.816798865471,0.0], [-20992.823456530015,9517.416843696026,0.0]])),
                         "QPQ5": survey_element(np.array([[-23329.778140908842,10844.180117790922,0.0], [-23287.611953749027,9295.383755414225,0.0]]))}

        self.nuin_to_bpd1 = self.nuin_to_bpd1 - self.startpoint
        self.bpd1_to_bpd2 = self.bpd1_to_bpd2 - self.startpoint
        self.bpd2_to_arc =  self.bpd2_to_arc -  self.startpoint
        for key, val in self.mag_objs.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec-self.startpoint
            self.mag_objs[key].name = key

        #drop the 3rd component
        self.nuin_to_bpd1 = self.nuin_to_bpd1[:, :2]
        self.bpd1_to_bpd2 = self.bpd1_to_bpd2[:, :2]
        self.bpd2_to_arc = self.bpd2_to_arc[:, :2]
        for key, val in self.mag_objs.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec[:, :2]

        print("after dropping 3rd component nuin to bpd1 ", self.nuin_to_bpd1)
       #rotate to point along (1, 0, 0)
        initial_vec = normvec(self.nuin_to_bpd1[1])
        print("initial vec ", initial_vec)
        rotangle = np.atan(initial_vec[1]/initial_vec[0])-np.pi
        print("rotating through angle ", rotangle)
        self.rotate(rotangle)

        self.svectors.append(self.nuin_to_bpd1[1] - self.nuin_to_bpd1[0])
        self.sstart.append(self.nuin_to_bpd1[0])

        self.svectors.append(self.bpd1_to_bpd2[1] - self.bpd1_to_bpd2[0])
        self.sstart.append(self.bpd1_to_bpd2[0])
        
        self.svectors.append(self.bpd2_to_arc[1] - self.bpd2_to_arc[0])
        self.sstart.append(self.bpd2_to_arc[0])
        print("svectors ", self.svectors)
        print("sstart ", self.sstart)

        self.get_intersections()
        self.get_survey()

    def draw_beamline(self):
        plt.plot(self.nuin_to_bpd1[:,0], self.nuin_to_bpd1[:,1])
        plt.plot(self.bpd1_to_bpd2[:,0], self.bpd1_to_bpd2[:,1])
        plt.plot(self.bpd2_to_arc[:,0], self.bpd2_to_arc[:,1])
        for key, val in self.mag_objs.items():
            plt.plot(val.nominal_normvec[:,0], val.nominal_normvec[:,1])

        plt.show()
    
    def rotate(self, theta):
        self.nuin_to_bpd1 = rotmat2d(theta).dot(self.nuin_to_bpd1.T).T
        self.bpd1_to_bpd2 = rotmat2d(theta).dot(self.bpd1_to_bpd2.T).T
        self.bpd2_to_arc = rotmat2d(theta).dot(self.bpd2_to_arc.T).T
        for key, val in self.mag_objs.items():
            self.mag_objs[key].nominal_normvec = rotmat2d(theta).dot(val.nominal_normvec.T).T
      
    def get_intersections(self):
        for key, val in self.mag_objs.items():
            #first check against nuin to bpd1
            summed_len = 0
            for i in range(len(self.svectors)):
                mu = get_intersection_len(self.sstart[i], self.svectors[i], val.nominal_normvec[0], val.nominal_normvec[1]-val.nominal_normvec[0])
                print(mu)
                if mu <= 1 and mu >=0:
                    val.s_nom = summed_len + mu*np.linalg.norm(self.svectors[i])
                    val.segment_id = i
                    val.segment = [np.append(self.sstart[i], 0.0), np.append(self.svectors[i], 0.0)]
                    val.nominal_center = self.sstart[i] + mu*self.svectors[i]
                    val.nominal_center = np.append(val.nominal_center, 0.0)
                    print("found intersection for ", key, " in segment ", i, " at length ", summed_len + mu*np.linalg.norm(self.svectors[i]))
                    break
                else:
                    summed_len += np.linalg.norm(self.svectors[i])
    def get_survey(self):
        self.line = strip_whitespace(pd.read_csv("../fujii-san.csv", header=0, skipinitialspace=True))
        self.survey = read_excel("03_2022_Neutrino.xlsx")
        print(self.survey)
        for key, value in self.mag_objs.items():

            self.mag_objs[key].nominal_offset = np.array(string_to_list(self.line.loc[self.line['element'] == key, 'survey_offset'].iloc[0]), dtype=np.float64)
#            print(f'survey offset for {key} is {self.mag_objs[key].nominal_offset}')
            for col in ['x2022', 'y2022', 'h2022']:
                self.mag_objs[key].survey_points.append(self.survey.loc[self.survey['name'] == key[1:]+'1', col].iloc[0])
            for col in ['x2022', 'y2022', 'h2022']:
                self.mag_objs[key].survey_points.append(self.survey.loc[self.survey['name'] == key[1:]+'2', col].iloc[0])

            self.mag_objs[key].survey_points = np.array(self.mag_objs[key].survey_points).reshape(2, 3)

#            print(f'survey points for {key} set to {self.mag_objs[key].survey_points}')

        #need to start to align the survey data to the nominal beamline, let's take BPV1 and QPQ2 (the most extreme magnets in the first straight stretch)
        bpv1_estimate = self.mag_objs['BPV1'].survey_points[0] - np.array([0., 0., self.mag_objs['BPV1'].nominal_offset[0,2]]) #subtract just the vertical component
        for mag in self.mag_objs.values():
            mag.survey_points = mag.survey_points - np.array([bpv1_estimate, bpv1_estimate])
#            print(f'changed {mag.name} survey points to {mag.survey_points}')

#        print(f'survey points for {key} set to {self.mag_objs[key].survey_points}')


        qpq2_estimate = self.mag_objs['QPQ2'].survey_points[0] - np.array([0., 0., self.mag_objs['QPQ2'].nominal_offset[0,2]]) #subtract just the vertical component
        bpv1_estimate = self.mag_objs['BPV1'].survey_points[0] - np.array([0., 0., self.mag_objs['BPV1'].nominal_offset[0,2]]) #subtract just the vertical component

        s_estimate = qpq2_estimate - bpv1_estimate
        #rotate in xy
        angle = np.pi+np.atan(s_estimate[1]/s_estimate[0]) #we know its in the bottom left quadrant
#        print(f'angle = {angle} for direction {s_estimate}')
        for mag in self.mag_objs.values():
            mag.rotate_survey('xy', -angle)
            mag.calculate_survey_center(mag.segment[1])
#            print(f'rotated survey data to {mag.survey_points}')

        #almost perfect at this point ~ 1mm precision over 50m of beamline
        #optional: repeat this rotation using the xz plane

        #now need to define a zero position, use BPV1 center from the survey and set it to the nominal position
        bpv1_survey_center = cp.deepcopy(self.mag_objs['BPV1'].survey_center)
        bpv1_nominal_center = self.mag_objs['BPV1'].nominal_center
        for mag in self.mag_objs.values():
            mag.shift_survey(bpv1_nominal_center - bpv1_survey_center)
            print(f'misalign for {mag.name} is {mag.survey_center - mag.nominal_center}')
        
        self.bline = beamline()
        print(self.bline.bpds[0].center, self.mag_objs['BPV1'].nominal_center)







#class to take the two survey points per element and determine the element center
class holder:
    point = []
    center = []
    offset = []
    set_center = False
    s = []
    name = []
    longitudinal = False
    def __init__(self, points, longitudinal=False, offsets=None, center=None, s=np.array([1., 0., 0.])):
        self.s = s.T/np.linalg.norm(s)
        self.longitudinal = longitudinal
        if(isinstance(points, pd.DataFrame)):
            self.point = np.array([points['x2022'].tolist(), points['y2022'].tolist(), points['h2022'].tolist()])
            self.point = [self.point[:,i] for i in range(self.point.shape[1])]
            print(points)
            self.name = points['name'].iloc[0]
        elif(isinstance(points[0], pd.DataFrame)):
            self.point = [np.array([pt['x2022'].iloc[0], pt['y2022'].iloc[0], pt['h2022'].iloc[0]]) for pt in points]
            self.name = points[0]['name'].iloc[0]
        elif(isinstance(points[0], list)):
            self.point = [cp.deepcopy(np.array(pt)) for pt in points]
            self.name = 'NA'
        elif(isinstance(points[0], np.ndarray)):
            self.point = [cp.deepcopy(pt) for pt in points]
            self.name = 'NA'
        elif(isinstance(points, holder)):
            self.point = cp.deepcopy(points.point)
            self.name = cp.deepcopy(points.name)
            self.center = cp.deepcopy(points.center)
        if center is not None:
            self.center = center
            self.set_center = True
        if offsets is not None:
            if(isinstance(offsets[0], np.ndarray)):
                self.offset = cp.deepcopy(offsets)
            elif(isinstance(offsets, list)):
                self.offset = [np.array(ele) for ele in offsets]
            elif(isinstance(offsets, str)):
                self.offset = np.array(string_to_list(offsets))
#                self.offset = [self.offset[:, i] for i in range(self.offset.shape[1])]
            else:
                print(offsets)
                raise Exception("Offset type undetermined ", type(offsets))
            self.calculate_center()
        else:
            self.offset = np.nan

    def __sub__(self, other):
        return holder([self.point[i] - other.point[i] for i in range(len(self.point))], center = self.get_center() - other.get_center())
    def __repr__(self):
        return f"class=holder, name='{self.name}', point='{self.point}', center='{self.center}'"

    def get_center(self):
       return self.center

    def binary_search(self, func, invtanth, R, s, conelen, l, r):
        precision = 1e-8
        val = 99999
        left = [l, s.dot(func(invtanth, R, conelen, l))]
        right = [r, s.dot(func(invtanth, R, conelen, r))]
        while(np.abs(val) > precision):
            testpoint = (left[0]+right[0])/2.
            val = s.dot(func(invtanth, R, conelen, testpoint))
            if(val<0):
                if(left[1]<0):
                    left = [testpoint, val]
                else:
                    right = [testpoint, val]
            elif(val>0):
                if(left[1]>0):
                    left = [testpoint, val]
                else:
                    right = [testpoint, val]
        if(np.abs(right[1]) < np.abs(left[1])):
            return right[0]
        else:
            return left[0]

    def calculate_center(self):
        #first get the vector between the two measured points
        vec = self.point[1] - self.point[0]
        offsetvec = self.offset[1] - self.offset[0]
        vecnorm = vec / np.linalg.norm(vec)
        offsetvecnorm = offsetvec / np.linalg.norm(offsetvec)
        #find a rotation matrix that takes offsetvec and rotates it to point along vec
        R = get_rotation_matrix(offsetvecnorm, vecnorm)
        p1c = R.dot(-self.offset[0]) #vector from point[0] to center
        cos_open_angle = vecnorm.dot(vec+p1c) / np.linalg.norm(vec+p1c)
        cone_opening_angle = np.acos(cos_open_angle)
        invtanopen = 1./np.tan(cone_opening_angle)
        conelen = np.linalg.norm(vec+p1c)
        p1cnorm = p1c/ np.linalg.norm(p1c)

        R2 = get_rotation_matrix(np.array([0., 0., 1.0]), -vec)

        locus = lambda invtan, R, conelen, ang: self.point[0] - (self.point[1] + conelen * R.dot(normvec([np.cos(ang), np.sin(ang), invtan]).T))

        vertical = np.array([0., 0., 1.])
        s = self.s
        if(self.longitudinal):
            s = np.cross(vertical, self.s)

        best_phi = self.binary_search(locus, invtanopen, R2, s, conelen, 0., np.pi)
        best_vector = self.point[0] + locus(invtanopen, R2, conelen, best_phi)
        if(best_vector[2]>self.point[0][2]): #got the solution above the reference points, but we want the one below
            best_phi = self.binary_search(locus, invtanopen, R2, s, conelen, np.pi, 2.*np.pi)
            best_vector = self.point[0] + locus(invtanopen, R2, conelen, best_phi)
            if(best_vector[2]>self.point[0][2]): #got the solution above the reference points, but we want the one below
                raise Exception("couldnt get a solution below the reference points")
        ## Create a figure
        #fig = plt.figure()
        ## Add a 3D subplot
        #ax = fig.add_subplot(111, projection='3d')
        #con = ax.scatter(center_points[:,0], center_points[:,1], center_points[:, 2])
        #reference1 = ax.scatter(self.point[0][0], self.point[0][1], self.point[0][2])
        #reference2 = ax.scatter(self.point[1][0], self.point[1][1], self.point[1][2])
        #best = ax.scatter(best_vector[0], best_vector[1], best_vector[2])
        #ax.set_xlabel('x(mm)')
        #ax.set_ylabel('y(mm)')
        #ax.set_zlabel('z(mm)')
        #ax.set_aspect('equal', 'box')
        #plt.show()
        #print()
        #print("Using vector ", R2.dot(normvec([np.cos(best_phi), np.sin(best_phi), invtanth])))
        self.center = best_vector
        self.set_center=True

    def project_along_beamline(self, vec):
        return self.point[0].dot(vec)

    def rotate(self, R):
        self.point = [R.dot(pt) for pt in self.point]
        self.center = R.dot(self.center)

    def set_zero(self, offset):
        self.point = [(pt + offset) for pt in self.point]
        self.center += offset
    def get_curvilinear_position(self):
        #first check if we're after BPDs
        after_bpds = False
        if(self.center[0]>end_of_bpd2):
            after_bpds = True
            start_s = end_of_bpd2

#class to store the survey data 
class beamline:
    line = []
    ssems=[]
    qpqs=[] 
    bpds =[]
    survey = []
    s = []
    def __init__(self):
        self.line = strip_whitespace(pd.read_csv("../fujii-san.csv", header=0, skipinitialspace=True))
        self.survey = read_excel("03_2022_Neutrino.xlsx")
        
        self.ssems, self.bpds, self.qpqs = parse_holder(self.survey, self.line, ssem1s)
        self.line['s_start'] = self.line['length'].shift().cumsum()
#        self.line['survey_offset_list'] = string_to_list(self.line['survey_offset'])
#        print(self.line['survey_offset_list'])
        self.map_surv_to_line()
        self.get_s()

    def map_surv_to_line(self):
        self.line['survey'] = np.nan
        for ss in self.ssems:
            self.line['survey'] = self.line.apply(lambda row: ss.center if row['element'] == ss.name[:-1] else row['survey'], axis=1)
        for qp in self.qpqs:
            self.line['survey'] = self.line.apply(lambda row: qp.center if row['element'][1:] == qp.name[:-1] else row['survey'], axis=1)
        for bp in self.bpds:
            self.line['survey'] = self.line.apply(lambda row: bp.center if row['element'][1:] == bp.name[:-1] else row['survey'], axis=1)

    def get_s(self):

        s = np.array(self.line.loc[self.line['element']=='SSEM5', 'survey']) - np.array(self.line.loc[self.line['element']=='SSEM4', 'survey'])
        datum = np.array(self.line.loc[self.line['element']=='SSEM4', 'survey'].iloc[0])
        datum[2] = 0. #just use ssem4 position for xy position, not height
        s[0][2] = 0.
        s = normvec(s[0])
        self.bendangle = np.atan(s[1]/s[0])
        bpd1idx = self.line.index[self.line['element'] == 'BPD1'][0]
        self.line['s_dir'] = self.line.apply(lambda row: s if row.name > bpd1idx else np.array([1., 0., 0.]), axis=1)
        self.line['angl'] = self.line.apply(lambda row: self.bendangle if row.name > bpd1idx else 0., axis=1)
        self.line['datum'] = self.line.apply(lambda row: datum if row.name > bpd1idx else np.array([0., 0., 0.]), axis=1)
       

        self.line['s_survey'] = self.line.apply(lambda row: np.array(row.s_dir).dot(np.array(row.survey)-np.array(row.datum)) if not np.isnan(row['survey']).any() else row['survey'], axis=1)

        self.line['misalign'] = self.line['survey']-self.line['datum'] - (self.line['s_survey']*self.line['s_dir'])
        #height shift is now correct but need some extra work for x and y after bending

       # self.line['shift'] = self.line.apply(lambda row: np.array(row.shift).dot(np.array(row.survey)-np.array(row.datum)) if not np.isnan(row['shift']).any() else np.array([0., 0., 0.]), axis=1)

        self.line['misalign'] = self.line.apply(lambda row: np.array([[np.cos(row.angl), np.sin(row.angl), 0.],[-np.sin(row.angl), np.cos(row.angl), 0.],[0., 0., 1.]]).dot(row.misalign) if not np.isnan(row.misalign).any() else np.nan, axis=1)

        #assume the difference between the x, y position of ssem4 and the s position of the start of BPD1 is the difference in s between ssem4 and bpd1, error ~0.3mm in s
        bpd1_s = self.line.loc[self.line['element'] == 'BPD1', 's_survey'].iloc[0]

        bpd1_ssem4_diff = np.array(self.line.loc[self.line['element'] == 'SSEM4', 'survey'].iloc[0] - [bpd1_s, 0., 0.])
        bpd1_ssem4_s_diff = np.linalg.norm(bpd1_ssem4_diff[:2])

        bpd1_bpd2_diff = np.array(self.line.loc[self.line['element'] == 'BPD2', 'survey'].iloc[0] - [bpd1_s, 0., 0.])
        bpd1_bpd2_s_diff = np.linalg.norm(bpd1_bpd2_diff[:2])
        self.line.loc[self.line['element'] == 'BPD2', 's_survey'] = bpd1_bpd2_s_diff + bpd1_s

        ssem4idx = self.line.index[self.line['element'] == 'SSEM4'][0]
        self.line['s_survey'] = self.line.apply(lambda row: row['s_survey'] + bpd1_ssem4_s_diff + bpd1_s if row.name >= ssem4idx else row['s_survey'], axis=1)        

        x = []
        y = []
        print("SSEM Position (mm): \t Nominal \t Survey \t Difference")
        for idx, row in self.line.iterrows():
            if not np.isnan(row['s_survey']):# and row['type'] == 'ssem':

                nominal_s_start = row['s_start']+float(row['mark'])-0.5*float(row['polelength'])
                x.append(nominal_s_start)
                y.append(nominal_s_start-row.s_survey)
                print(row['element'],',\t\t\t',"{:.3f}".format(nominal_s_start), ',\t',"{:.3f}".format(row['s_survey']), ',\t',"{:.3f}".format(nominal_s_start-row['s_survey']))
#        plt.scatter(x, y)
#        plt.show()
        self.line = self.line.drop(['survey', 's_dir', 'angl', 'datum'], axis=1)



class BeamlinePrinter:
    def __init__(self, line, kv, filename, primaries_only=False):
        self.beamline = line
        self.kvals = kv
        self.file = open(filename, "w")
        self.s = 0
        self.blmID = 1
        self.primaries_only=primaries_only
        self.line = []

    ######################################################
    #beam properties
    ######################################################

    def print_beam_sad(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0=-0.0004542861867420988*m,
	  Xp0=2.5505774863429934e-05,
      emitx=0.084015*mm*mrad,
      betx=37.098*m,
      alfx=-2.41877,
      dispx=0.423734*m,
      dispxp=0.0719639,

      Y0=-0.00023643145922826628*m,
	  Yp0=7.751587096049882e-05,
      emity=0.0695782*mm*mrad,
      bety=5.45*m,
      alfy=0.178,
      dispy=0.0*m,
      dispyp=0.,

      kineticEnergy=30*GeV;\n\n''')

    def print_beam(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0=-0.0004542861867420988*m,
      Xp0=2.5505774863429934e-05,
      emitx=0.084015*mm*mrad,
      betx=39.8775555502868*m,
      alfx=-0.568424785315,
      dispx=0.423734*m,
      dispxp=0.0719639,

      Y0=-0.00023643145922826628*m,
      Yp0=7.751587096049882e-05,
      emity=0.0695782*mm*mrad,
      bety=6.4519061397745*m,
      alfy=0.35934594606,
      dispy=0.0*m,
      dispyp=0.,

      kineticEnergy=30*GeV;\n\n''')

    def print_beam_0910580(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0=0.0*m,
      Xp0=0.0,
      emitx=0.19*mm*mrad,
      betx=35.835*m,
      alfx=-2.3704,
      dispx=0.443*m,
      dispxp=0.074,

      Y0=0.0*m,
      Yp0=0.0,
      emity=0.157*mm*mrad,
      bety=7.369*m,
      alfy=0.064,
      dispy=0.0*m,
      dispyp=0.,

      kineticEnergy=30*GeV;\n\n''')

    def print_beam_0910216(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0=0.0*m,
      Xp0=0.0,
      emitx=0.075116*mm*mrad,
      betx=37.75916*m,
      alfx=-2.33673231,
      dispx=0.033185*m,
      dispxp=0.00337915,

      Y0=0.0*m,
      Yp0=0.0,
      emity=0.06006*mm*mrad,
      bety=5.5537*m,
      alfy=0.19780927,
      dispy=0.0*m,
      dispyp=0.,

      kineticEnergy=30*GeV;\n\n''')

    def print_beam_sadfit_0910216(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0=-0.5*mm,
      Xp0=3.5e-5,
      emitx=0.0610768*mm*mrad,
      betx=37.098*m,
      alfx=-2.4187,
      dispx=0.42373*m,
      dispxp=0.07196,
                        
      Y0=-0.2*mm,
      Yp0=7.8e-5,
      emity=0.05976*mm*mrad,
      bety=5.45*m,
      alfy=0.178,
      dispy=0.0*m,
      dispyp=0.,

      kineticEnergy=30*GeV;\n\n''')
        
    def print_beam_from_file(self, filename):
        self.file.write('''beam, particle="proton",
      kineticEnergy=30*GeV,
      distrType="bdsimsampler:entry",''')
        self.file.write('distrFile="'+filename+'";\n')

    def print_halo(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="halo",

      X0=-0.0004542861867420988*m,
      Xp0=2.5505774863429934e-05,
      emitx=0.084015*mm*mrad,
      betx=39.8775555502868*m,
      alfx=-0.568424785315,

      Y0=-0.00023643145922826628*m,
      Yp0=7.751587096049882e-05,
      emity=0.0695782*mm*mrad,
      bety=6.4519061397745*m,
      alfy=0.35934594606,
      haloNSigmaXInner=4.0,
      haloNSigmaYInner=4.0,                       
      haloNSigmaXOuter=8.0,
      haloNSigmaYOuter=8.0,
      haloPSWeightFunction="flat",

      kineticEnergy=30*GeV;\n\n
        ''')


    def endl(self):
        self.file.write(";\n")


    ######################################################
    #now the physical elements
    ######################################################

    def print_aperture(self, row):
        self.file.write(', apertureType="'+str(row.aperture_type)+'"')
        if(row.aperture_type == 'circular'): 
            self.file.write(', aper1=' + str(0.5*row.aperture_x) + '*mm')
        elif(row.aperture_type == 'rectangular'):
            self.file.write(', aper1=' + str(0.5*row.aperture_x) + '*mm, aper2=' + str(0.5*row.aperture_y) + '*mm')

    def print_drift(self, row, name, driftlen):
        self.line.append(name)
        self.file.write(name + ': drift, l=' + str(driftlen)+ '*mm')
        self.print_aperture(row)
        self.print_xsec_bias('vacuum')

    def print_bend_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt)+', B='+str(self.kvals[row.element])+'*T')
        if(misalignments):
            self.file.write(', offsetX='+str(row.misalign[1])+', offsetY='+str(row.misalign[2])+', offsetZ='+str(row.misalign[0]))
        if(not geometry):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('vacuum')
        self.endl()

    def print_quad_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, tilt='+str(row.tilt)+', k1='+str(self.kvals[row.element]))
        if(misalignments):
            self.file.write(', offsetX='+str(row.misalign[1])+', offsetY='+str(row.misalign[2])+', offsetZ='+str(row.misalign[0]))
        if(not geometry):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('vacuum')
        self.endl()

    def print_ssem(self, row, thickness):
        first_driftlen = row.mark
        name = row.element + "_udrift"
        self.print_drift(row, name, row.mark)
        self.endl()
        misalign = row.misalign
        if(np.isnan(row.misalign).any()):
            misalign = [0., 0., 0.]
        self.print_target(row.element, str(thickness)+"*mm", "G4_Ti", row.aperture_x, misalign)
        name = row.element + "_ddrift"
        self.print_drift(row, name, row.length - float(row.mark) - thickness) 
        self.endl()

    def print_target(self, name, length, material, hWidth, misalign):
        self.line.append(name)
        if(misalignments):
            self.file.write(name+": target, l="+str(length)+", material=\""+material+"\", horizontalWidth="+str(hWidth)+"*mm, offsetX="+str(misalign[1])+"*mm, offsetY="+str(misalign[2])+"*mm, offsetZ="+str(misalign[0])+"*mm")
        else:
            self.file.write(name+": target, l="+str(length)+", material=\""+material+"\", horizontalWidth="+str(hWidth)+"*mm")
        self.print_xsec_bias('material')
        self.endl()

    def print_dump(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': dump, horizontalWidth='+str(row.aperture_x)+'*mm, l='+str(row.length)+'*mm')
        self.endl()

    def print_blm(self, reference, dx, dy, ds, orientation):
        self.file.write('blm_'+reference+'_'+str(self.blmID)+': blm, scoreQuantity="chrg eDep", geometryType="cylinder", blm1=100*mm, blm2=30*mm, blmMaterial="Al",')
        self.file.write('referenceElement="'+reference+'", x='+str(dx)+'*mm, y='+str(dy)+'*mm, s='+str(ds)+'*mm')
        if(orientation=='perp'):
            self.file.write(', theta=1.570796, psi=1.570796')
        #TODO maybe this also needs xsec biasing?
        self.endl()
        self.blmID+=1

    def print_tunnel(self):
        self.file.write('''
option, buildTunnel = 1,
tunnelType="rectangular",
tunnelOffsetX = 100*cm,
tunnelOffsetY = 50*cm,
tunnelAper1 = 150*cm,
tunnelAper2 = 150*cm,
tunnelThickness = 30*cm,
buildTunnelFloor = 0,
tunnelSoilThickness = 2*m;\n\n''')


    ######################################################
    #now for the more abstract things
    ######################################################
    def print_xsec_bias(self, typ):
        if(bias_physics):
            if(typ == 'vacuum'):
                self.file.write(', bias="vacBias"')
            elif(typ == 'material'):
                self.file.write(', biasMaterial="matBias"')

    def print_field(self, row, ndim):
         self.file.write(row.element+'field: field, type="bmap'+str(ndim)+'d", bScaling=1.0, magneticFile="bdsim'+str(ndim)+'d:../magnet_responses/'+row.element+'.dat"')
         self.endl()
       
    def print_fieldmapgeom(self, row, ndim):
        self.line.append(row.element)
        self.print_field(row, ndim)
        self.file.write(row.element+': element, geometryFile="gdml:../'+row.element+'.gdml", fieldAll="'+row.element+'field", l='+str(row.length)+'*mm')
        self.print_xsec_bias('vacuum')
        self.endl()

    def print_fieldmap(self, row, magtype, ndim):
        self.line.append(row.element)
        self.print_field(row, ndim)
        self.file.write(row.element+': '+magtype+', fieldVacuum="'+row.element+'field", l='+str(row.length)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt))
        if(not geometry):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('vacuum')
        self.endl()

    def print_blms(self, row):
        dx = string_to_list(row.blm_offset_x)
        dy = string_to_list(row.blm_offset_y)
        ds = [str(float(entry)+row.polelength/2.0) for entry in string_to_list(row.blm_offset_s)]  #TODO THIS ONLY WORKS FOR PURE QUAD/DIPOLE FIELDMAP WILL BREAK THIS
        orientation = string_to_list(row.blm_orientation)
        for i in range(len(dx)):
            self.print_blm(row.element, dx[i], dy[i], ds[i], orientation[i])

    def print_physics(self, physics_list):
        self.file.write('option, physicsList =' + '"' + physics_list + '";\n')

    def sample_ssems(self):
        for element in self.line:
            if(re.match('SSEM[0-9]$', element)):
                self.file.write('sample, range='+element+';\n')

    ######################################################
    #main function to decide which functions to call based on element type
    ######################################################

    def print(self):
        self.file.write('chrg: scorer, type="cellcharge";\n')
        self.file.write('eDep: scorer, type="depositeddose";\n')
        previously_drift = False
        prev_row = []
        driftlen = 0.0
        for row in self.beamline.itertuples():
            
            if(row.type != 'drift' and previously_drift and merge_drifts):
                #we have reached a non-drift element, so flush the drift components to the file
                self.print_drift(prev_row, prev_row.element, driftlen)
                self.endl()
                driftlen = 0.0

            if(row.type in ['rbend', 'sbend', 'quadrupole']):
                self.file.write('\n') #give some breathing room
                self.print_drift(row, row.element+'_driftu', float(row.mark)-row.polelength/2.)
                self.endl()
                if(row.type == 'quadrupole'):
                    self.print_quad_magnet(row)
                else:
                    self.print_bend_magnet(row)
                self.print_drift(row, row.element+'_driftd', row.length - (float(row.mark)+row.polelength/2.))
                self.endl()
                self.file.write('\n') #give some breathing room

            elif(re.match('fieldmap.*', row.type)):
                ndim = extract_number(row.type)
                if(re.match('.*geom', row.type)):
                    self.print_fieldmapgeom(row, ndim)
                else:
                    magtypes = ['sbend', 'rbend', 'quadrupole']
                    for magtype in magtypes:
                        print(row.type)
                        if(re.match('.*'+magtype, row.type)):
                            print("found match "+magtype)
                            self.print_fieldmap(row, magtype, ndim)

#            elif(row.type == 'fieldmap3dquad'):
#                self.print_fieldmap(row, 'drift') #TODO left these three here in case we want to add geometry dependant on magtype
#            elif(row.type == 'fieldmap3drbend'):
#                self.print_fieldmap(row, 'drift')
#            elif(row.type == 'fieldmap3dsbend'):
#                self.print_fieldmap(row, 'drift')

            elif(row.type == 'ssem'):
                self.print_ssem(row, 15e-3)
            elif(row.type == 'wsem'):
                self.print_ssem(row, 1e-4)
            elif(row.type == 'dump'):
                self.print_dump(row)
            else: #otherwise, i.e. drift or collimator
                if(not merge_drifts):
                    self.print_drift(row, row.element, row.length)
                    self.endl()
                else:
                    driftlen += row.length


            prev_row = row
            if(row.type == 'drift'):
                previously_drift = True

            self.s += row.length
            if(row.blm and enable_blms):
                self.print_blms(row)
            self.file.write("! s=" + str(self.s) + "\n")

        #print the beamline elements in a line
        self.file.write('\nl0: line = (')
        if(self.primaries_only):
            self.line = [self.line[0]]
        for idx, element in enumerate(self.line):
            if idx == len(self.line)-1:
                self.file.write(element +');\n')
            else:
                self.file.write(element +',\n')

        self.file.write('\nuse, period=l0;\n')
        #self.print_beam()
        if(beam_from_file):
            self.print_beam_from_file("../run_sadfit_0910216_10k.root")
        elif(beam_halo):
            self.print_halo()
        else:
            self.print_beam_sadfit_0910216()

        if(print_tunnel):
            self.print_tunnel()
        if(print_physics):
            self.print_physics('g4FTFP_BERT')
        self.file.write('option, nturns=1;\n')
        if(sample_all):
            self.file.write('sample, all;\n')
        elif(sample_ssem):
            self.sample_ssems()
        if(sample_entry):
            self.file.write('sample, range=entry;\n')
        if(print_vacuum):
            self.file.write('option, vacuumPressure='+str(vacuum_pressure)+';\n')
        if(bias_physics): #if we're biasing the physics cross-section flag=2 means only applies to primaries
            self.file.write('vacBias: xsecBias, particle="proton", proc="all", xsecfact=1, flag=2;\n')
            self.file.write('matBias: xsecBias, particle="proton", proc="all", xsecfact=1, flag=2;\n')


def get_holder(surv, basename, points, longitudinal, offsets, s=None):
    names = [basename+str(pt) for pt in points]
    if(s is None):
        return holder(surv[surv['name'].isin(names)], longitudinal, offsets)
    else:
        return holder(surv[surv['name'].isin(names)], longitudinal, offsets, s)


def read_excel(filename):
    surv = pd.ExcelFile(filename)
    surv = surv.parse('基準点成果表 (側壁fit) ', skiprows=3)
    surv.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters
    surv.columns = ['index', 'ID', 'name', 'x2022', 'y2022', 'h2022', 'x2017', 'y2017', 'h2017', 'x2014', 'y2014', 'h2014', 'dx', 'dy', 'dh', 'comment']
    surv = surv.drop(['x2014', 'y2014', 'h2014'], axis=1)
    surv = surv.drop(['x2017', 'y2017', 'h2017'], axis=1)
    surv = surv.drop(['dx', 'dy', 'dh'], axis=1)
    return surv

def rotate_surveyxy(surv, theta):
    tmp = cp.deepcopy(surv)
    surv['x2022'] = np.cos(theta) * tmp['x2022'] - np.sin(theta) * tmp['y2022']
    surv['y2022'] = np.sin(theta) * tmp['x2022'] + np.cos(theta) * tmp['y2022']
    return surv

def rotate_surveyxz(surv, theta):
    tmp = cp.deepcopy(surv)
    R = np.array([[np.cos(theta), np.sin(theta)],
                [-np.sin(theta), np.cos(theta)]])
    surv['x2022'] = R[0,0] * tmp['x2022'] + R[0,1] * tmp['h2022']
    surv['h2022'] = R[1,0] * tmp['x2022'] + R[1,1] * tmp['h2022']
    return surv

def parse_holder(surv, line, ssem1s):

#    ssem_offset = line[re.match('ssem[0-9]*', line['element'])] ##SSEM offsets are the same for all SSEMs (for now, not true in FF)
    bpd_offset = line[line['element'].str.contains(r'^BP', regex=True)]
    bpd_names = bpd_offset['element']
    bpd_offset = bpd_offset.set_index('element')['survey_offset'].to_dict()

    qpq_offset = line[line['element'].str.contains(r'^QP', regex=True)]

    print("qpq_offset_type ", type(qpq_offset['survey_offset']))

    qpq_names = qpq_offset['element']
    qpq_offset = qpq_offset.set_index('element')['survey_offset'].to_dict()
    ssem_offset = np.array([[0., 0., 225.5],[0., 330., 225.5]]) #offset in x,y,h where x is along the beamline,y is horisontal and h is vertical
#    bpd_offset =  np.array([[0., 0., 500.],[3000., 0., 500.]]) #offset in x,y,h where x is along the beamline,y is horisontal and h is vertical
#    qpq_offset =  np.array([[0., 0., 500.],[3000., 0., 500.]]) #offset in x,y,h where x is along the beamline,y is horisontal and h is vertical

    #due to the issue of multiple solutions, described in the rotation matrix we must make an assumption:
    #components are placed such that the vertical axis of their vertical face is perpendicular to the beamline 
    #essentially the points are on top of the components, so allow them to hang from these points and be pulled vertical by gravity


    #SSEM1/2 are God
    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    ssem2 = [surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']]
   
    #estimate initial beamline direction taking ssem2[0]-ssem1[0]
    s = np.array([ssem2[0][i].iloc[0] - ssem1[0][i].iloc[0] for i in ['x2022', 'y2022', 'h2022']])

    #get the xy rotation
    theta = np.atan(s[1]/s[0])
    theta = -np.pi-theta #we know that the survey is in the bottom-right quadrant
    #rotate the whole survey to point along 1, 0, z
    surv = rotate_surveyxy(surv, theta)

    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    ssem2 = [surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']]
    s = np.array([ssem2[0][i].iloc[0] - ssem1[0][i].iloc[0] for i in ['x2022', 'y2022', 'h2022']])

    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    ssem1_center = holder(ssem1, False, ssem_offset, s=s).get_center()


    #coords relative to ssem1[0]
    surv.loc[:,'x2022'] = surv['x2022']-ssem1_center[0]
    surv.loc[:,'y2022'] = surv['y2022']-ssem1_center[1]
    surv.loc[:,'h2022'] = surv['h2022']-ssem1_center[2]

    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    ssem2 = [surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']]
    s = np.array([ssem2[0][i].iloc[0] - ssem1[0][i].iloc[0] for i in ['x2022', 'y2022', 'h2022']])

    ssem1_center = holder([surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']], False, ssem_offset, s=s).get_center()
    ssem2_center = holder([surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']], False, ssem_offset, s=s).get_center()
    #again rotate to point along 1., 0., z but this time use the centers
    theta = np.atan((ssem2_center[1] - ssem1_center[1])/(ssem2_center[0]-ssem1_center[0]))
    surv = rotate_surveyxy(surv, -theta)


    #and again in the xz plane
    ssem1_center = holder([surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']], False, ssem_offset, s=s).get_center()
    ssem2_center = holder([surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']], False, ssem_offset, s=s).get_center()
    s = ssem2_center - ssem1_center

    theta = np.atan((ssem2_center[2] - ssem1_center[2])/(ssem2_center[0]-ssem1_center[0]))
    surv = rotate_surveyxz(surv, theta)



    ssem1_center = holder([surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']], False, ssem_offset, s=s).get_center()
    ssem2_center = holder([surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']], False, ssem_offset, s=s).get_center()
    s = ssem2_center - ssem1_center

    ang = bpd1_angle+bpd2_angle
    R = np.array([[np.cos(ang), np.sin(ang), 0.],
                [-np.sin(ang), np.cos(ang), 0.],
                [0., 0., 1.]])

    s_after_bpd = R.dot(s)
    ssem_svec = [s if i<=3 else s_after_bpd for i in range(1,10)]
    qpq_svec = {'QPQ1': s,
             'QPQ2': s,
             'QPQ3': s_after_bpd,
             'QPQ4': s_after_bpd,
             'QPQ5': s_after_bpd}
    
    bpd_svec = {'BPV1': s,
             'BPH2': s,
             'BPD1': s,             #CERN TODO these should be something between s and s_after_bpd
             'BPD2': s_after_bpd,   #CERN TODO 
             'BPV2': s_after_bpd,
             'BPH3': s_after_bpd}
   


    #qpq_svec = [s if i<=2 else s_after_bpd for i in range(1,6)]

    ssems = [get_holder(surv, 'SSEM'+str(id), [1, 2], False, ssem_offset, s=ssem_svec[id-1]) for id in range(1,10)]
    qpqs = [get_holder(surv, id[1:], [1, 2], True, qpq_offset[id], s=qpq_svec[id]) for id in qpq_names]
    bpds = [get_holder(surv, id[1:], [1, 2], True, bpd_offset[id], s=bpd_svec[id]) for id in bpd_names]

    #NOTE SSEMs now 0 indexed!!!

    beamline_dir = ssems[1].get_center() - ssems[0].get_center()
    beamline_dir /= np.linalg.norm(beamline_dir)
#    print("beamline dir ", beamline_dir)

    offset = np.array([ssem1s, 0., 0.])
    [ss.set_zero(offset) for ss in ssems]
    [bp.set_zero(offset) for bp in bpds]
    [qp.set_zero(offset) for qp in qpqs]
    return ssems, bpds, qpqs

def plot_points(pts, lab, xcomp, ycomp, zcomp=None, ax=None):
    x = []
    y = []
    z = []
    for i in range(len(pts)):
        for j in range(len(pts[i].point)):
            x.append(pts[i].point[j][xcomp])
            y.append(pts[i].point[j][ycomp])
            if(zcomp is not None):
                z.append(pts[i].point[j][zcomp])
    
    if(ax is not None):
        ax.scatter(x, y, z, label=lab)
    else:
        plt.scatter(x, y, label=lab)

def plot_centers(hld, lab, xcomp, ycomp, zcomp=None, ax=None):
    x = []
    y = []
    z = []
    for i in range(len(hld)):
        x.append(hld[i].get_center()[xcomp])
        y.append(hld[i].get_center()[ycomp])
        if(zcomp is not None):
            z.append(hld[i].get_center()[zcomp])
    if(ax is not None):
        ax.scatter(x, y, z, label=lab, marker='x')
    else:
        plt.scatter(x, y, label=lab, marker='x')
      

if __name__ == '__main__':

    nom = nominalBeamline()
    nom.draw_beamline()
    exit(1)

    beamline_with_misalign = beamline()
    bline = beamline_with_misalign.line

    #run 910216
    vec_magset = [0 ,
    -15 ,
    520 ,
    0 ,
    485 ,
    1140 ,
    1191 ,
    408 ,
    15 ,
    354 ,
    -13 ,
    423 ]
    
    
    magnet_response = strip_whitespace(pd.read_csv("../kicurve.csv", header=0, skipinitialspace=True))
    
    
    #copy of the magnet mapping in SAD
    magset = {}
    magset["BPV1"] = vec_magset[0]
    magset["BPH2"] = vec_magset[1]
    magset["QPQ1"] = vec_magset[2]
    magset["QPQ2"] = vec_magset[4]
    magset["BPD1"] = vec_magset[5]
    magset["BPD2"] = vec_magset[6]
    magset["QPQ3"] = vec_magset[7]
    magset["BPV2"] = vec_magset[8]
    magset["QPQ4"] = vec_magset[9]
    magset["BPH3"] = vec_magset[10]
    magset["QPQ5"] = vec_magset[11]
    
    
    kvals = {}
    
    for magnet in magset:
        mag_df = magnet_response[magnet_response['element'] == magnet]
        kvals[magnet] = np.interp(magset[magnet], mag_df['current'], mag_df['kval'])
        zero_field = np.interp(0, mag_df['current'], mag_df['kval'])
        if magnet[0] == 'B': #bending magnets
            kvals[magnet] = -(kvals[magnet]-zero_field) * (proton_momentum/0.2998)  / (0.001*bline.loc[bline['element'] == magnet].iloc[0]['polelength'])
        else:
            kvals[magnet] = (kvals[magnet]-zero_field) / (0.001*bline.loc[bline['element'] == magnet].iloc[0]['polelength'])
        if abs(kvals[magnet]) < 1e-8: #if the strength is zero bdsim will treat it as a drift so force it to be non-zero
            kvals[magnet] = 1e-8
    
    
    #kvals['BPD1'] = -1.15329
    #kvals['BPD2'] = -1.14018
    #kvals['QPQ4'] = -0.0518735 #set the QPQ4 val to that in the fake fieldmap
    
    mag_df = magnet_response[magnet_response['element'] == 'BPH3']
    #plt.scatter(mag_df['current'], mag_df['kval'])
    #plt.show()
    
    
    #df_tmp = magnet_response[magnet_response['element'] == row.element]
    
    
    prnt = BeamlinePrinter(bline, kvals, "test.gmad", primaries_only=generate_primaries)
    prnt.print()


#    exit(1)

###################just for plotting purposes###################
    survey = read_excel("03_2022_Neutrino.xlsx")
    ssems, bpds, qpqs = parse_holder(survey, bline, ssem1s)


    s = np.linspace(0, 50000, 10000)
    vec_curvilinear_coords = np.vectorize(curvilinear_coords)
    x, y = vec_curvilinear_coords(s)

    mu = np.linspace(0, 50000, 1000)


    plot_points(ssems, 'SSEM Survey', 0, 2)
    plot_points(bpds, 'BPD Survey', 0, 2)
    plot_points(qpqs, 'QPQ Survey', 0, 2)
    plot_centers(ssems, 'SSEM center', 0, 2)
    plot_centers(bpds, 'BPD center', 0, 2)
    plot_centers(qpqs, 'QPQ center', 0, 2)
    plt.xlabel('x(mm)')
    plt.ylabel('z(mm)')
    plt.legend()
    plt.show()


    plot_points(ssems, 'SSEM', 0, 1)
    plot_points(bpds, 'BPD', 0, 1)
    plot_points(qpqs, 'QPQ', 0, 1)
    plot_centers(ssems, 'SSEM center', 0, 1)
    plot_centers(bpds, 'BPD center', 0, 1)
    plot_centers(qpqs, 'QPQ center', 0, 1)
    plt.xlabel('x(mm)')
    plt.ylabel('y(mm)')
    plt.legend()
    plt.show()


    fig = plt.figure()
    # Add a 3D subplot
    ax = fig.add_subplot(111, projection='3d')
    plot_points(ssems, 'SSEM', 0, 1, 2, ax)
    plot_centers(ssems, 'SSEM center', 0, 1, 2, ax)
    plot_points(bpds, 'BPD', 0, 1, 2, ax)
    plot_centers(bpds, 'BPD center', 0, 1, 2, ax)
    plot_points(qpqs, 'QPQ', 0, 1, 2, ax)
    plot_centers(qpqs, 'QPQ center', 0, 1, 2, ax)
    ax.set_xlabel('x(mm)')
    ax.set_ylabel('y(mm)')
    ax.set_zlabel('z(mm)')

    plt.legend()
    #ax.set_aspect('equal', 'box')
    plt.show()




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

##beam loss configuration
#print_tunnel=False
#print_physics=False
#sample_all=True
#sample_ssem=True
#sample_entry=False  #####WARNING MUST BE FALSE WHEN FITTING OTHERWISE ENTRY WILL BE TREATED AS SSEM1!!!
#beam_from_file = False
#beam_halo = False
#enable_blms = True
#geometry = True
#misalignments=True
#bias_physics=False
#print_vacuum=False
#merge_drifts=False


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
    
def string_to_list(s):
    s = s.replace(' ', ',')
    s = s.replace('][', '],[')
    elements = ast.literal_eval(s)
    return elements

def string_to_string_list(s):
    return s.strip("[]").split()

def read_excel(filename):
    surv = pd.ExcelFile(filename)
    surv = surv.parse('基準点成果表 (側壁fit) ', skiprows=3)
    surv.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters
    surv.columns = ['index', 'ID', 'name', 'x2022', 'y2022', 'h2022', 'x2017', 'y2017', 'h2017', 'x2014', 'y2014', 'h2014', 'dx', 'dy', 'dh', 'comment']
    surv = surv.drop(['x2014', 'y2014', 'h2014'], axis=1)
    surv = surv.drop(['x2017', 'y2017', 'h2017'], axis=1)
    surv = surv.drop(['dx', 'dy', 'dh'], axis=1)
    return surv




class survey_element:
    def __init__(self, nom_normvec):
        self.nominal_normvec = cp.deepcopy(nom_normvec)
        self.s_nom = 0.
        self.segment = []
        self.segment_id = -1.
        self.nominal_center = []
        self.survey_points = []
        self.survey_center = []
        self.survey_diff = []
        self.nominal_offset = []
        self.name = []
        self.magnum = -1
        self.segment_start_s = []
        self.angle_diff = 0.
        self.angle = 0.

    def rotate_survey(self, axis, theta):
        c = np.cos(theta)
        s = np.sin(theta)
        if(axis == 'xy'):
            R = np.array([[c, -s, 0.],[s, c, 0.], [0., 0., 1.]])
        elif(axis == 'xz'):
            R = np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])
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
        precision = 1e-10
        val = 99999
        left = [l, func(l, R, horisontal)]
        right = [r, func(r, R, horisontal)]
        niter = 0
        while(np.abs(val) > precision and niter < 1000):
            niter += 1
            testpoint = (left[0]+right[0])/2.
            val = func(testpoint, R, horisontal)
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

    def calculate_survey_center(self):

        survey_point_diff = cp.deepcopy(self.survey_points[1] - self.survey_points[0])
        nominal_surv_point_diff = cp.deepcopy(-self.nominal_offset[1] + self.nominal_offset[0]) #these numbers are center to point vector so sort of backwards from the survey

        R = get_rotation_matrix(nominal_surv_point_diff, survey_point_diff) #rotation to take the points from the idealised space to the survey space

        p1c = R.dot(-self.nominal_offset[0]) #direction to go from survey point 1 to the center
        #there is an extra DOF, rotation about the line connecting the two survey points traces out a cone of valid p1c vectors
        #how far along the magnet direction does the p1c line go?
        len_along_magdir = p1c.dot(normvec(-survey_point_diff))
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
        #get resulting center point
        self.survey_center = self.survey_points[0] + len_along_magdir*normvec(survey_point_diff) + circle_radius*R2.dot(np.array([np.cos(best_angle), np.sin(best_angle), 0.0]))
        #alternate method as a cross-check
        center = cp.deepcopy(self.survey_points[0])
        center -= self.nominal_offset[0][0]*normvec(survey_point_diff) #move along vector connecting survey points by amount nominal diff
        center -= self.nominal_offset[0][1]*horisontal
        center -= self.nominal_offset[0][2]*vertical

        if(np.linalg.norm(center - self.survey_center) > 2e-1): ##if the two methods differ by more than 0.2mm (tested to be below threshold normally)
            print(f'WARNING centers calculated through two methods differ above threshold for {self.name}, {center} vs {self.survey_center}')
            print(f'Difference is {center-self.survey_center} be concerned if the second or third component is the driver of the difference')


class nominalBeamline:
    svectors = []
    sstart = []
    slen = []
    bline = []
    offsets = []
    mag_objs = []
    def __init__(self): #a bunch of vectors extraced from the .dwg files provided by Fujii-san
        self.line = strip_whitespace(pd.read_csv("../fujii-san.csv", header=0, skipinitialspace=True)) #the csv containing beampipe properties
        self.line['s_start'] = self.line['length'].shift().cumsum()
        for element in self.line.itertuples():
#            print(element)
            if(re.match('SSEM[0-9]', element.element)):
                print(element.element, element.s_start+float(element.mark))



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
        ctr = 0
        for key, val in self.mag_objs.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec-self.startpoint
            self.mag_objs[key].name = key
            self.mag_objs[key].magnum = ctr
            ctr += 1

        #drop the 3rd component
        self.nuin_to_bpd1 = self.nuin_to_bpd1[:, :2]
        self.bpd1_to_bpd2 = self.bpd1_to_bpd2[:, :2]
        self.bpd2_to_arc = self.bpd2_to_arc[:, :2]
        for key, val in self.mag_objs.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec[:, :2]

       #rotate to point along (1, 0, 0)
        initial_vec = normvec(self.nuin_to_bpd1[1])
        rotangle = np.atan(initial_vec[1]/initial_vec[0])-np.pi
        self.rotate(rotangle)

        self.svectors.append(self.nuin_to_bpd1[1] - self.nuin_to_bpd1[0])
        self.sstart.append(self.nuin_to_bpd1[0])

        self.svectors.append(self.bpd1_to_bpd2[1] - self.bpd1_to_bpd2[0])
        self.sstart.append(self.bpd1_to_bpd2[0])
        
        self.svectors.append(self.bpd2_to_arc[1] - self.bpd2_to_arc[0])
        self.sstart.append(self.bpd2_to_arc[0])

        self.get_intersections()
        self.get_survey()
        self.get_offsets()
        self.match_misalignments()


    def draw_beamline_s(self, axis):
        vars = ['s (m)', 'Horisontal (mm)', 'Vertical (mm)']
        nom_s = []
        for key, val in self.mag_objs.items():
            plt.scatter((val.s_nom+val.offsets[0])/1000., val.offsets[axis], label=val.name+' survey')

            nom_s.append(val.s_nom/1000.)
        plt.plot(nom_s, np.zeros(len(nom_s)), label='Nominal')
        plt.xlabel(vars[0])
        plt.xlim(0, 1.05*np.max(nom_s))
        plt.ylabel(vars[axis])
        plt.legend(loc='upper right')
        plt.show()

       

    def draw_beamline(self, axis1, axis2):
        labels = ['x (mm)', 'y (mm)', 'z(mm)']

        nuin_to_bpd1_3d = np.hstack((self.nuin_to_bpd1, np.zeros((self.nuin_to_bpd1.shape[0], 1))))
        bpd1_to_bpd2_3d = np.hstack((self.bpd1_to_bpd2, np.zeros((self.bpd1_to_bpd2.shape[0], 1))))
        bpd2_to_arc_3d = np.hstack((self.bpd2_to_arc, np.zeros((self.bpd2_to_arc.shape[0], 1))))

        plt.plot(nuin_to_bpd1_3d[:,axis1], nuin_to_bpd1_3d[:,axis2])
        plt.plot(bpd1_to_bpd2_3d[:,axis1], bpd1_to_bpd2_3d[:,axis2])
        plt.plot(bpd2_to_arc_3d[:,axis1], bpd2_to_arc_3d[:,axis2])
        for key, val in self.mag_objs.items():
            nominal_normvec_3d = np.hstack((val.nominal_normvec, np.zeros((val.nominal_normvec.shape[0], 1))))

#            plt.plot(nominal_normvec_3d[:,axis1], nominal_normvec_3d[:,axis2], label=val.name)
            plt.scatter(val.survey_points[:,axis1], val.survey_points[:,axis2], label=val.name)
            plt.scatter(val.survey_center[axis1], val.survey_center[axis2], label=val.name)
        plt.legend(loc='upper right')
        plt.xlabel(labels[axis1])
        plt.ylabel(labels[axis2])
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
                if mu <= 1 and mu >=0:
                    val.segment_start_s = summed_len
                    val.s_nom = summed_len + mu*np.linalg.norm(self.svectors[i])
                    val.segment_id = i
                    val.segment = [np.append(self.sstart[i], 0.0), np.append(self.svectors[i], 0.0)]
                    val.nominal_center = self.sstart[i] + mu*self.svectors[i]
                    val.nominal_center = np.append(val.nominal_center, 0.0)
                    if(mu*np.linalg.norm(self.svectors[i]) < 1.0): #bending magnet (within 1mm of the change point)
                        val.angle = np.acos(normvec(self.svectors[i]).dot(normvec(self.svectors[i-1]))) * np.sign(normvec(self.svectors[i-1])[1] - normvec(self.svectors[i])[1]) # negative if y increases
                    elif(mu*np.linalg.norm(self.svectors[i]) > np.linalg.norm(self.svectors[i]) - 1.0): #bending magnet (within 1mm of the change point)
                        val.angle = np.acos(normvec(self.svectors[i]).dot(normvec(self.svectors[i+1]))) * np.sign(normvec(self.svectors[i])[1] - normvec(self.svectors[i+1])[1]) # negative if y increases
                    print("found intersection for ", key, " in segment ", i, " at length ", summed_len + mu*np.linalg.norm(self.svectors[i]), " with angle ", val.angle)
                    break
                else:
                    summed_len += np.linalg.norm(self.svectors[i])
    def get_survey(self):
        self.line = strip_whitespace(pd.read_csv("../fujii-san.csv", header=0, skipinitialspace=True))
        self.survey = read_excel("03_2022_Neutrino.xlsx")
        for key, value in self.mag_objs.items():

            self.mag_objs[key].nominal_offset = np.array(string_to_list(self.line.loc[self.line['element'] == key, 'survey_offset'].iloc[0]), dtype=np.float64)
            for col in ['x2022', 'y2022', 'h2022']:
                self.mag_objs[key].survey_points.append(self.survey.loc[self.survey['name'] == key[1:]+'1', col].iloc[0])
            for col in ['x2022', 'y2022', 'h2022']:
                self.mag_objs[key].survey_points.append(self.survey.loc[self.survey['name'] == key[1:]+'2', col].iloc[0])

            self.mag_objs[key].survey_points = np.array(self.mag_objs[key].survey_points).reshape(2, 3)


    def get_offsets(self):

        #need to start to align the survey data to the nominal beamline, let's take BPV1 and QPQ2 (the most extreme magnets in the first straight stretch)
        bpv1_estimate = self.mag_objs['BPV1'].survey_points[0] - np.array([0., 0., self.mag_objs['BPV1'].nominal_offset[0,2]]) #subtract just the vertical component
        for mag in self.mag_objs.values():
            mag.survey_points = mag.survey_points - np.array([bpv1_estimate, bpv1_estimate])

        qpq2_estimate = self.mag_objs['QPQ2'].survey_points[0] - np.array([0., 0., self.mag_objs['QPQ2'].nominal_offset[0,2]]) #subtract just the vertical component
        bpv1_estimate = self.mag_objs['BPV1'].survey_points[0] - np.array([0., 0., self.mag_objs['BPV1'].nominal_offset[0,2]]) #subtract just the vertical component

        s_estimate = qpq2_estimate - bpv1_estimate
        #approximate rotatation in xy, technically not necessary but makes understanding/debug so much easier
        angle = np.pi+np.atan(s_estimate[1]/s_estimate[0]) #we know its in the bottom left quadrant
        for mag in self.mag_objs.values():
            mag.rotate_survey('xy', -angle)
            mag.calculate_survey_center()

        #almost perfect at this point ~ 2-3mm precision over 50m of beamline
        #repeat this rotation using the xz plane and using the survey centers rather than the survey points themselves
        #use largest lever arm possible, in the vertical this is BPV1 and QPQ5
        angle = np.atan((self.mag_objs['QPQ5'].survey_center[2]-self.mag_objs['BPV1'].survey_center[2])/(self.mag_objs['QPQ5'].survey_center[0]-self.mag_objs['BPV1'].survey_center[0]))
        for mag in self.mag_objs.values():
            mag.rotate_survey('xz', -angle)
            mag.calculate_survey_center()
        #and again for xy in y there are bending magnets BPD1,2 so can only use the line to QPQ2 for the zero direction
        angle = np.atan((self.mag_objs['QPQ2'].survey_center[1]-self.mag_objs['BPV1'].survey_center[1])/(self.mag_objs['QPQ2'].survey_center[0]-self.mag_objs['BPV1'].survey_center[0]))
        for mag in self.mag_objs.values():
            mag.rotate_survey('xy', -angle)
            mag.calculate_survey_center()

        #sadly BPD1 doesn't bend by the correct amount relative to the survey, while we could technically offset everything its better to just adjust BPD2 bending magnitude
        dr = self.mag_objs['QPQ5'].survey_center - self.mag_objs['QPQ3'].survey_center 
        survey_angle = np.atan(dr[1]/dr[0])
        nominal_angle = np.atan(self.mag_objs['QPQ3'].segment[1][1] / self.mag_objs['QPQ3'].segment[1][0])
        self.mag_objs['BPD2'].angle += nominal_angle - survey_angle
        #now change the segment angle do not change the segment start position 

        for key, val in self.mag_objs.items():
            if(val.magnum > self.mag_objs['BPD2'].magnum): #after BPD2
                self.mag_objs[key].segment[1] = np.linalg.norm(val.segment[1]) * normvec(np.array([np.cos(survey_angle), np.sin(survey_angle), 0.0]))
                #need to recalculate nominal center
                self.mag_objs[key].nominal_center = val.segment[0] + (val.s_nom - val.segment_start_s) * normvec(val.segment[1])

        #now need to define a zero position, use BPV1 center from the survey and set it to the nominal position
        bpv1_survey_center = cp.deepcopy(self.mag_objs['BPV1'].survey_center)
        bpv1_nominal_center = cp.deepcopy(self.mag_objs['BPV1'].nominal_center)
        vertical = np.array([0., 0., 1.])
        for mag in self.mag_objs.values():
            mag.shift_survey(bpv1_nominal_center - bpv1_survey_center)
            #now need misalignment in beamline coords
            #along beamline
            offset = mag.survey_center - mag.nominal_center
            along_segment = (mag.survey_center - mag.segment[0]).dot(normvec(mag.segment[1])) 
            offset_s = along_segment - mag.s_nom + mag.segment_start_s 
            vertical_offset = offset.dot(vertical)
            horisontal_vector = -normvec(np.cross(mag.segment[1], vertical))
            horisontal_offset = horisontal_vector.dot(offset)

            mag.offsets = [offset_s, horisontal_offset, vertical_offset]
        
    def match_misalignments(self):
        self.line['misalign'] = np.nan

        for key, value in self.mag_objs.items():
            self.line['misalign'] = self.line.apply(lambda row: value.offsets if row['element'] == key else row['misalign'], axis=1)
            if(misalignments):
                self.line['angle'] = self.line.apply(lambda row: value.angle if row['element'] == key else row['angle'], axis=1)


    def print_beamline(self, kv, filename):
        self.printer = BeamlinePrinter(self.line, kv, filename)
        self.printer.print()


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

    def print_aperture(self, row, vertical=False):
        self.file.write(', apertureType="'+str(row.aperture_type)+'"')
        if(row.aperture_type == 'circular'): 
            self.file.write(', aper1=' + str(0.5*row.aperture_x) + '*mm')
        elif(row.aperture_type == 'rectangular'):
            if(vertical):
                self.file.write(', aper1=' + str(0.5*row.aperture_y) + '*mm, aper2=' + str(0.5*row.aperture_x) + '*mm')
            else:
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
            self.file.write(', offsetX='+str(row.misalign[1])+'*mm, offsetY='+str(row.misalign[2])+'*mm')
        if(not geometry):
            self.file.write(', magnetGeometryType="none"')
        vertical = False
        if(np.abs(row.tilt-np.pi/2.0) <0.01):
            vertical=True
        self.print_aperture(row, vertical)
        self.print_xsec_bias('vacuum')
        self.endl()

    def print_quad_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, tilt='+str(row.tilt)+', k1='+str(self.kvals[row.element]))
        if(misalignments):
            self.file.write(', offsetX='+str(row.misalign[1])+'*mm, offsetY='+str(row.misalign[2])+'*mm')
        if(not geometry):
            self.file.write(', magnetGeometryType="none"')
        vertical = False
        if(np.abs(row.tilt-np.pi/2.0) <0.01):
            vertical=True
        self.print_aperture(row, vertical)
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
            self.file.write(name+": target, l="+str(length)+", material=\""+material+"\", horizontalWidth="+str(hWidth)+"*mm, offsetX="+str(misalign[1])+"*mm, offsetY="+str(misalign[2])+"*mm")
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
        orientation = string_to_string_list(row.blm_orientation)
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
                if(misalignments):
                    self.print_drift(row, row.element+'_driftu', float(row.mark)-row.polelength/2.+row.misalign[0])
                else:
                    self.print_drift(row, row.element+'_driftu', float(row.mark)-row.polelength/2.)

                self.endl()
                if(row.type == 'quadrupole'):
                    self.print_quad_magnet(row)
                else:
                    self.print_bend_magnet(row)
                if(misalignments):
                    self.print_drift(row, row.element+'_driftd', row.length - (float(row.mark)+row.polelength/2.)-row.misalign[0])
                else:
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



if __name__ == '__main__':

    nom = nominalBeamline()

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
            kvals[magnet] = -(kvals[magnet]-zero_field) * (proton_momentum/0.2998)  / (0.001*nom.line.loc[nom.line['element'] == magnet].iloc[0]['polelength'])
        else:
            kvals[magnet] = (kvals[magnet]-zero_field) / (0.001*nom.line.loc[nom.line['element'] == magnet].iloc[0]['polelength'])
        if abs(kvals[magnet]) < 1e-3: #if the strength is zero bdsim will treat it as a drift so force it to be non-zero, if its too small the integrator will fall over however
            kvals[magnet] = 1e-3
    
    
    #kvals['BPD1'] = -1.15329
    
    
    nom.print_beamline(kvals, "test.gmad")


    nom.draw_beamline_s(1)
    nom.draw_beamline_s(2)
    nom.draw_beamline(0, 1)
    nom.draw_beamline(0, 2)

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import re
import copy as cp
import sys
import ast
import ROOT
from array import array 
#sys.path.append('../')
#import create_beamline
#matplotlib.use('qtagg')
pd.set_option('display.max_rows', 500000)
#pd.set_option('display.max_columns', 500000)
false=False
true=True

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
if(True):
    print_tunnel=False
    print_physics=False
    sample_all=False
    sample_ssem=True
    sample_entry=False  #####WARNING MUST BE FALSE WHEN FITTING OTHERWISE ENTRY WILL BE TREATED AS SSEM1!!!
    beam_from_file = False
    beam_halo = False
    use_previous_best_fit = False
    enable_blms = False
    geometry = False
    misalignments=True
    bias_physics=False
    print_vacuum=False
    merge_drifts=True

##beam orbit plot configuration
if(False):
    print_tunnel=False
    print_physics=False
    sample_all=True
    sample_ssem=False
    sample_entry=False  #####WARNING MUST BE FALSE WHEN FITTING OTHERWISE ENTRY WILL BE TREATED AS SSEM1!!!
    beam_from_file = False
    beam_halo = False
    use_previous_best_fit = True
    enable_blms = False
    geometry = False
    misalignments=True
    bias_physics=False
    print_vacuum=False
    merge_drifts=False

##beam loss configuration
if(False):
    print_tunnel=False
    print_physics=True
    sample_all=True
    sample_ssem=False
    sample_entry=False  #####WARNING MUST BE FALSE WHEN FITTING OTHERWISE ENTRY WILL BE TREATED AS SSEM1!!!
    beam_from_file = False
    beam_halo = False
    use_previous_best_fit = True
    enable_blms = True
    geometry = True
    misalignments=True
    bias_physics=True
    print_vacuum=True
    merge_drifts=False

#generate primaries configuration
if(False):
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
def rotmatxy3d(theta):
    return np.array([[np.cos(theta), np.sin(theta), 0],
                    [-np.sin(theta), np.cos(theta), 0],
                    [0., 0., 1.]])
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

    def rotate_nominal(self, theta):
        self.nominal_normvec = rotmat2d(theta).dot(self.nominal_normvec.T).T
        if isinstance(self.segment, np.ndarray): #don't do anything if segment hasnt been set
            
            print(f"rotating segment as well {self.segment}")
            self.segment = rotmatxy3d(theta).dot(self.segment.T).T
            self.calculate_nominal_center()
            print(f"rotated segment to {self.segment}")


    def shift_survey(self, change):
        self.survey_center += change
        self.survey_points += change

    def shift_nominal(self, change):
        print(f"adding {change} to {self.name}")
        self.nominal_normvec += change[:2] #CERN TODO only able to shift nominal normvec in the plane
        self.nominal_center += change
        self.segment[0] += change

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

    def calculate_nominal_center(self):
        self.nominal_center = self.segment[0] + (self.s_nom - self.segment_start_s) * normvec(self.segment[1])



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



        self.nuin_to_bpd1 = np.array([[27711.320296225927,10246.557442498312], [8826.078822922398,10998.629426709593]])
        self.startpoint_ps = self.nuin_to_bpd1[0,:]

        self.bpd1_to_bpd2 = np.array([[8826.078822922398,10998.629426709593], [4625.740813728273,11025.057932632815]])
        self.bpd2_to_arc = np.array([[4625.740813728273,11025.057932632815],[-26953.86986321407,10165.298827699111]])


        self.mag_objs_ps = {"BPV1": survey_element(np.array([[24832.38910462714,9576.830981585355,0.0], [24885.262495664367,10904.532021424293,0.0]])),
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


        ff_height_offset = -309.2167344454938 #the height in the drawing of QFQ1 center
        self.arc_to_fh1 = np.array([[17088.847808979302,8354.310785989503], [9691.14780897929,8364.639930467885]])
        self.startpoint_ff = self.arc_to_fh1[0,:]
        self.fh1_to_fh2 = np.array([[9691.14780897929,8364.639930467885], [-18598.159457284823,8364.639930467885]])
        self.fvd1_axis = np.array([[-8196.453655041289,8364.639930467885,-306.1337343460509-ff_height_offset],[-10570.220438481108,8364.639930467885,-336.2819644230731-ff_height_offset]])
        self.fq4_axis = np.array([[-11262.47885696281,8364.639930467885,-355.3463325089069-ff_height_offset],[-13161.86822774537,8364.639930467885,-403.51289931818155-ff_height_offset]])
        self.fvd2_axis = np.array([[-13920.76196973481,8364.639930467885,-417.21062514883306-ff_height_offset],[-15739.152191020708,8364.639930467885,-493.74168506447313-ff_height_offset]])

        print(f'nominal BFVD1 bend angle = {np.arccos(normvec(self.fq4_axis[0]-self.fq4_axis[1]).dot([1.0, 0.0, 0.0]))}')
        print(f'nominal BFVD2 bend angle = {2.0*np.arccos(normvec(self.fvd2_axis[0]-self.fvd2_axis[1]).dot(normvec(self.fq4_axis[0]-self.fq4_axis[1])))}')


        self.mag_objs_ff={"QFQ1": survey_element(np.array([[14491.769195548608,8807.02230689272,0.0], [14490.532107179635,7921.023170543049,0.0]])),
                         "BFV1": survey_element(np.array([[11790.204692362477,7351.24453761384,0.0], [11791.911689156808,9185.063540965684,0.0]])),
                         "BFH1": survey_element(np.array([[9691.14780897929,9185.063540965684,0.0], [9691.14780897929,7354.310785989502,0.0]])),
                         "BFV2": survey_element(np.array([[3730.8478089792916,9128.096944042545,0.0], [3730.8478089792916,7651.590156712228,0.0]])),
                         "QFQ2": survey_element(np.array([[780.8478089792907,9203.346478591457,0.0], [780.8478089792907,7615.581943857097,0.0]])),
                         "QFQ3": survey_element(np.array([[-3969.1521910207093,9203.346478591457,0.0], [-3969.1521910207093,7615.581943857097,0.0]])),
                         "BFH2": survey_element(np.array([[-7019.15219101996,9803.376606143274,0.0], [-7019.152191021339,7136.859410107129,0.0]])),
#                         "BFVD1": survey_element(np.array([[-9430.288218548048,9194.639930556094, 0.0], [-9430.288218548048,7594.639930467868,0.0]])), ##CERN TODO removing these so the nominal angle is applied
                         "QFQ4": survey_element(np.array([[-11473.32213158561-715.0650297,8320.339930467722,0.0], [-11473.32213158561-715.0650297,8389.939930467724,0.0]]))}#, #some maffs needed for x
#                         "BFVD2": survey_element(np.array([[-15709.697351297618,7594.639930503842,0.0], [-15709.697351297618,9194.639930592066,0.0]]))}

        self.mag_objs = {**self.mag_objs_ps, **self.mag_objs_ff}


        self.nuin_to_bpd1 = self.nuin_to_bpd1 - self.startpoint_ps
        self.bpd1_to_bpd2 = self.bpd1_to_bpd2 - self.startpoint_ps
        self.bpd2_to_arc =  self.bpd2_to_arc -  self.startpoint_ps

        self.fvd1_axis -= np.append(self.arc_to_fh1[0,:], 0.0) 
        self.fq4_axis -= np.append(self.arc_to_fh1[0,:], 0.0) 
        self.fvd2_axis -= np.append(self.arc_to_fh1[0,:], 0.0)

        self.arc_to_fh1 = self.arc_to_fh1 - self.startpoint_ff
        self.fh1_to_fh2 = self.fh1_to_fh2 - self.startpoint_ff

        #lets get the intersections for the FF bending down section while all the vectors still have the same y component
        #apologies for the horrible element accesses but necessary to choose the x, z axes
        self.fh2_to_fvd1 = get_intersection_len(np.append(self.fh1_to_fh2[0,0], 0.0), np.append(self.fh1_to_fh2[1,0]-self.fh1_to_fh2[0,0], 0.0), self.fvd1_axis[0][[0,2]], (self.fvd1_axis[1]-self.fvd1_axis[0])[[0,2]])
        print(f"fh2 to fvd1 is at length {self.fh2_to_fvd1} along fh1 to fh2 line {self.fh1_to_fh2[1]-self.fh1_to_fh2[0]}")
        print(f"point position = {self.fh1_to_fh2[0] + self.fh2_to_fvd1*(self.fh1_to_fh2[1]-self.fh1_to_fh2[0])}")
        print(f"for reference fh1_to_fh2 {self.fh1_to_fh2}")
        self.fh2_to_fvd1 = np.array([np.append(self.fh1_to_fh2[0] + self.fh2_to_fvd1*(self.fh1_to_fh2[1]-self.fh1_to_fh2[0]), 0.0), self.fvd1_axis[1]])

        print(f"fh2 to fvd1 is given by {self.fh2_to_fvd1}")
        self.fvd1_to_fvd2 = get_intersection_len(self.fh2_to_fvd1[0][[0,2]], self.fh2_to_fvd1[1][[0,2]] - self.fh2_to_fvd1[0][[0,2]], self.fq4_axis[0][[0,2]], self.fq4_axis[1][[0,2]] - self.fq4_axis[0][[0,2]])
        print(f"length along fq4_axis is {self.fvd1_to_fvd2}")

        ctr = 0

        self.startpoint_ps = np.append(self.startpoint_ps, 0.0)
        for key, val in self.mag_objs_ps.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec-self.startpoint_ps
            self.mag_objs[key].name = key
            self.mag_objs[key].magnum = ctr
            ctr += 1
        
        self.startpoint_ff = np.append(self.startpoint_ff, 0.0)
        for key, val in self.mag_objs_ff.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec-self.startpoint_ff
            self.mag_objs[key].name = key
            self.mag_objs[key].magnum = ctr
            ctr += 1
           

        #drop the 3rd component
        self.nuin_to_bpd1 = self.nuin_to_bpd1[:, :2]
        self.bpd1_to_bpd2 = self.bpd1_to_bpd2[:, :2]
        self.bpd2_to_arc = self.bpd2_to_arc[:, :2]



        for key, val in self.mag_objs.items():
            self.mag_objs[key].nominal_normvec = val.nominal_normvec[:, :2]

        #rotate prep section to point along (1, 0, 0)
        initial_vec = normvec(self.nuin_to_bpd1[1])
        rotangle = np.arctan(initial_vec[1]/initial_vec[0])-np.pi
        self.nuin_to_bpd1 = rotmat2d(rotangle).dot(self.nuin_to_bpd1.T).T
        self.bpd1_to_bpd2 = rotmat2d(rotangle).dot(self.bpd1_to_bpd2.T).T
        self.bpd2_to_arc = rotmat2d(rotangle).dot(self.bpd2_to_arc.T).T
        self.rotate_nominal(rotangle, self.mag_objs_ps)
        
        #rotate ff section to point along (0, 1, 0)
        initial_vec = normvec(self.arc_to_fh1[1])
        rotangle = np.arctan(initial_vec[1]/initial_vec[0])+0.5*np.pi

        self.arc_to_fh1 = rotmat2d(rotangle).dot(self.arc_to_fh1.T).T
        self.fh1_to_fh2 = rotmat2d(rotangle).dot(self.fh1_to_fh2.T).T
        self.fvd1_axis = rotmatxy3d(rotangle).dot(self.fvd1_axis.T).T
        self.fq4_axis = rotmatxy3d(rotangle).dot(self.fq4_axis.T).T
        self.fvd2_axis = rotmatxy3d(rotangle).dot(self.fvd2_axis.T).T
        self.rotate_nominal(rotangle, self.mag_objs_ff)
        [print(self.mag_objs[key].nominal_normvec) for key, val in self.mag_objs_ff.items()]

#        get_intersection(self.fh1_to_fh2, self.fvd1_axis)



        self.svectors.append(self.nuin_to_bpd1[1] - self.nuin_to_bpd1[0])
        self.sstart.append(self.nuin_to_bpd1[0])

        self.svectors.append(self.bpd1_to_bpd2[1] - self.bpd1_to_bpd2[0])
        self.sstart.append(self.bpd1_to_bpd2[0])
        
        self.svectors.append(self.bpd2_to_arc[1] - self.bpd2_to_arc[0])
        self.sstart.append(self.bpd2_to_arc[0])

        self.svectors.append(self.arc_to_fh1[1] - self.arc_to_fh1[0])
        self.sstart.append(self.arc_to_fh1[0])

        self.svectors.append(self.fh1_to_fh2[1] - self.fh1_to_fh2[0])
        self.sstart.append(self.fh1_to_fh2[0])

#        self.sstart.append(fh1_to_fh2_intersection_bfvd1)



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
            if(val.name == 'QFQ1'):
                break
            nominal_normvec_3d = np.hstack((val.nominal_normvec, np.zeros((val.nominal_normvec.shape[0], 1))))

            plt.quiver(val.segment[0][axis1], val.segment[0][axis2], val.segment[1][axis1], val.segment[1][axis2], angles='xy', scale_units='xy', scale=1, width=0.004)
            plt.plot(nominal_normvec_3d[:,axis1], nominal_normvec_3d[:,axis2], label=val.name)
            plt.scatter(val.survey_points[:,axis1], val.survey_points[:,axis2], label=val.name)
            plt.scatter(val.survey_center[axis1], val.survey_center[axis2], label=val.name)
        plt.legend(loc='upper right')
        plt.xlabel(labels[axis1], fontsize=13)
        plt.ylabel(labels[axis2], fontsize=13)
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        plt.show()
    
    def rotate_nominal(self, theta, iterlist):
        for key, val in iterlist.items():
            self.mag_objs[key].rotate_nominal(theta)

    def shift_nominal(self, shift, iterlist):
        for key, val in iterlist.items():
            self.mag_objs[key].shift_nominal(shift)

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
                    val.segment = np.array([np.append(self.sstart[i], 0.0), np.append(self.svectors[i], 0.0)])
                    val.nominal_center = self.sstart[i] + mu*self.svectors[i]
                    val.nominal_center = np.append(val.nominal_center, 0.0)
                    if(mu*np.linalg.norm(self.svectors[i]) < 1.0): #bending magnet (within 1mm of the change point)
                        val.angle = np.arccos(normvec(self.svectors[i]).dot(normvec(self.svectors[i-1]))) * np.sign(normvec(self.svectors[i-1])[1] - normvec(self.svectors[i])[1]) # negative if y increases
                    elif(mu*np.linalg.norm(self.svectors[i]) > np.linalg.norm(self.svectors[i]) - 1.0): #bending magnet (within 1mm of the change point)
                        val.angle = np.arccos(normvec(self.svectors[i]).dot(normvec(self.svectors[i+1]))) * np.sign(normvec(self.svectors[i])[1] - normvec(self.svectors[i+1])[1]) # negative if y increases
                    print("found intersection for ", key, " in segment ", i, " at length ", summed_len + mu*np.linalg.norm(self.svectors[i]), " with angle ", val.angle, " mu ", mu)
                    print(val.segment)
                    break
                else:
                    summed_len += np.linalg.norm(self.svectors[i])

    def get_survey(self):
        self.line_ps = strip_whitespace(pd.read_csv("../fujii-san.csv", header=0, skipinitialspace=True))
        self.line_arc = strip_whitespace(pd.read_csv("../fujii-san_arc.csv", header=0, skipinitialspace=True))
        self.line_ff = strip_whitespace(pd.read_csv("../fujii-san_FF.csv", header=0, skipinitialspace=True))
        self.line = pd.concat([self.line_ps, self.line_arc, self.line_ff], ignore_index=True)
        self.survey = read_excel("03_2022_Neutrino.xlsx")
        for key, value in self.mag_objs.items():
            self.mag_objs[key].nominal_offset = np.array(string_to_list(self.line.loc[self.line['element'] == key, 'survey_offset'].iloc[0]), dtype=np.float64)
            for col in ['x2022', 'y2022', 'h2022']:
                self.mag_objs[key].survey_points.append(self.survey.loc[self.survey['name'] == key[1:]+'1', col].iloc[0])
            for col in ['x2022', 'y2022', 'h2022']:
                self.mag_objs[key].survey_points.append(self.survey.loc[self.survey['name'] == key[1:]+'2', col].iloc[0])

            self.mag_objs[key].survey_points = np.array(self.mag_objs[key].survey_points).reshape(2, 3)

    def get_offsets(self):

        #need to start to align the survey data to the nominal beamline, let's take QPQ1 and QPQ2 (the most extreme quads in the first straight stretch)
        QPQ1_estimate = self.mag_objs['QPQ1'].survey_points[0] - np.array([0., 0., self.mag_objs['QPQ1'].nominal_offset[0,2]]) #subtract just the vertical component
        for mag in self.mag_objs.values():
            mag.survey_points = mag.survey_points - np.array([QPQ1_estimate, QPQ1_estimate])

        qpq2_estimate = self.mag_objs['QPQ2'].survey_points[0] - np.array([0., 0., self.mag_objs['QPQ2'].nominal_offset[0,2]]) #subtract just the vertical component
        QPQ1_estimate = self.mag_objs['QPQ1'].survey_points[0] - np.array([0., 0., self.mag_objs['QPQ1'].nominal_offset[0,2]]) #subtract just the vertical component

        s_estimate = qpq2_estimate - QPQ1_estimate
        #approximate rotatation in xy, technically not necessary but makes understanding/debug so much easier
        angle = np.pi+np.arctan(s_estimate[1]/s_estimate[0]) #we know its in the bottom left quadrant
        for mag in self.mag_objs.values():
            mag.rotate_survey('xy', -angle)
            mag.calculate_survey_center()

        #almost perfect at this point ~ 2-3mm precision over 50m of beamline
        #repeat this rotation using the xz plane and using the survey centers rather than the survey points themselves
        #use largest lever arm possible, in the vertical this is QPQ1 and QPQ5
        angle = np.arctan((self.mag_objs['QPQ5'].survey_center[2]-self.mag_objs['QPQ1'].survey_center[2])/(self.mag_objs['QPQ5'].survey_center[0]-self.mag_objs['QPQ1'].survey_center[0]))
        for mag in self.mag_objs.values():
            mag.rotate_survey('xz', -angle)
            mag.calculate_survey_center()
        #and again for xy in y there are bending magnets BPD1,2 so can only use the line to QPQ2 for the zero direction
        angle = np.arctan((self.mag_objs['QPQ2'].survey_center[1]-self.mag_objs['QPQ1'].survey_center[1])/(self.mag_objs['QPQ2'].survey_center[0]-self.mag_objs['QPQ1'].survey_center[0]))
        for mag in self.mag_objs.values():
            mag.rotate_survey('xy', -angle)
            mag.calculate_survey_center()

        #sadly BPD1 doesn't bend by the correct amount relative to the survey, while we could technically offset everything its better to just adjust BPD2 bending magnitude
        dr = self.mag_objs['QPQ5'].survey_center - self.mag_objs['QPQ3'].survey_center 
        survey_angle = np.arctan(dr[1]/dr[0])
        nominal_angle = np.arctan(self.mag_objs['QPQ3'].segment[1][1] / self.mag_objs['QPQ3'].segment[1][0])
        self.mag_objs['BPD2'].angle += nominal_angle - survey_angle
        #now change the segment angle do not change the segment start position 

        for key, val in self.mag_objs.items():
            if(val.magnum > self.mag_objs['BPD2'].magnum and val.magnum<self.mag_objs['QFQ1'].magnum): #after BPD2 and before FF
                self.mag_objs[key].segment[1] = np.linalg.norm(val.segment[1]) * normvec(np.array([np.cos(survey_angle), np.sin(survey_angle), 0.0]))
                #need to recalculate nominal center
                self.mag_objs[key].calculate_nominal_center()# = val.segment[0] + (val.s_nom - val.segment_start_s) * normvec(val.segment[1])

        #even worse, we need to create a fake arc section, use the input vector connecting QFQ1 and BFH1
        ff_input = self.mag_objs['BFH1'].survey_center - self.mag_objs['QFQ1'].survey_center
        start_vector = self.mag_objs['QPQ2'].survey_center - self.mag_objs['QPQ1'].survey_center
        print("survey FF input vector", ff_input)
        print("survey PS input vector ", start_vector)
        ff_input_anglexy = np.arctan(ff_input[1]/ff_input[0])
        print("survey FF input angle ", ff_input_anglexy)

        #rotate the nominal normvecs
        self.rotate_nominal(0.5*np.pi-ff_input_anglexy, self.mag_objs_ff)
        #lets use the ps nominal positions as our reference frame, need to shift and rotate ff nominals into this system
        #already QFQ1 center is defined at 0, 0, 0 but need to rotate it to point along the line connecting QFQ1 to BFH1
        print(f"QFQ1 survey center {self.mag_objs['QFQ1'].survey_center}")
#        for key, val in self.mag_objs.items():
#            if(val.magnum >= self.mag_objs['QFQ1'].magnum): #after QFQ1
                #take QFQ1 survey position as the nominal position and add this position to all magnets after QFQ1
#                self.mag_objs[key].segment[0] += self.mag_objs['QFQ1'].survey_center
                #need to recalculate nominal center
#                self.mag_objs[key].calculate_nominal_center()# = val.segment[0] + (val.s_nom - val.segment_start_s) * normvec(val.segment[1])
       

        #now need to define a zero position, use QPQ1 center from the survey and set it to the nominal position
        qpq1_survey_center = cp.deepcopy(self.mag_objs['QPQ1'].survey_center)
        qpq1_nominal_center = cp.deepcopy(self.mag_objs['QPQ1'].nominal_center)
        for mag in self.mag_objs.values():
            mag.shift_survey(qpq1_nominal_center - qpq1_survey_center)

            
        self.shift_nominal(self.mag_objs['QFQ1'].survey_center - self.mag_objs['QFQ1'].nominal_center, self.mag_objs_ff)
        arc_qfq1_offset = (self.mag_objs['QFQ1'].s_nom - self.mag_objs['QFQ1'].segment_start_s)*normvec(self.mag_objs['QFQ1'].segment[1])
        print(f"arc to FQ1 vector {arc_qfq1_offset}")
        print(f"QFQ1 segment direction {self.mag_objs['QFQ1'].segment[1]}")
        print(f"QFQ1 to BFH1 vector {self.mag_objs['BFH1'].survey_center - self.mag_objs['QFQ1'].survey_center}")
#        for key, value in self.mag_objs_ff.items():
#                self.mag_objs[key].shift_nominal(self.mag_objs['QFQ1'].survey_center) ##TODO needs to be shifted by the difference between QFQ1 and FFinput point
        
        vertical = np.array([0., 0., 1.])
        for mag in self.mag_objs.values():
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
#            if(misalignments): #CERN TODO
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
        self.terminal_element = 'd1'

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
      emitx=0.075116e-6*m*rad,
      betx=37.75916*m,
      alfx=-2.33673231,
      dispx=0.033185*m,
      dispxp=0.00337915,

      Y0=0.0*m,
      Yp0=0.0,
      emity=0.06006e-6*m*rad,
      bety=5.5537*m,
      alfy=0.19780927,
      dispy=0.0*m,
      dispyp=0.,

      !sigmaE=0.3e-2,

      kineticEnergy=30*GeV;\n\n''')

    def print_beam_sadfit_0910216(self):
        self.file.write('''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0=-0.0005*m,
      Xp0=3.5e-5,
      emitx=0.0610768e-6*m*rad,
      betx=37.098*m,
      alfx=-2.4187,
      dispx=0.42373*m,
      dispxp=0.07196,
                        
      Y0=-0.0002*m,
      Yp0=7.8e-5,
      emity=0.05976e-6*m*rad,
      bety=5.45*m,
      alfy=0.178,
      dispy=0.0*m,
      dispyp=0.,

      !sigmaE=0.3e-2,
                        
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

    def print_beam(self):

      self.file.write(f'''\n\nbeam, particle="proton",
      distrType="gausstwiss",

      X0={kvals["X0"]}*m,
      Xp0={kvals["Xp0"]},
      emitx={kvals["emitx"]}*m*rad,
      betx={kvals["betx"]}*m,
      alfx={kvals["alfx"]},
      dispx=0.423734*m,
      dispxp=0.0719639,

      Y0={kvals["Y0"]}*m,
      Yp0={kvals["Yp0"]},
      emity={kvals["emity"]}*m*rad,
      bety={kvals["bety"]}*m,
      alfy={kvals["alfy"]},
      dispy=0.0*m,
      dispyp=0.,

      !sigmaE=0.3e-2,

      kineticEnergy=30*GeV;\n\n''')
      
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
        self.print_xsec_bias('')

    def print_bend_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt)+', B='+str(self.kvals[row.element])+'*T')
        if(misalignments):
            self.file.write(', offsetX='+str(row.misalign[1])+'*mm, offsetY='+str(row.misalign[2])+'*mm')
        if(not geometry or row.geometry=='none'):
            self.file.write(', magnetGeometryType="none"')
        vertical = False
        if(np.abs(row.tilt-np.pi/2.0) <0.01):
            vertical=True
        self.print_aperture(row, vertical)
        self.print_xsec_bias('')
        self.endl()

    def print_quad_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, tilt='+str(row.tilt)+', k1='+str(self.kvals[row.element]))
        if(misalignments):
            self.file.write(', offsetX='+str(row.misalign[1])+'*mm, offsetY='+str(row.misalign[2])+'*mm')
        if(not geometry or row.geometry == 'none'):
            self.file.write(', magnetGeometryType="none"')
        vertical = False
        if(np.abs(row.tilt-np.pi/2.0) <0.01):
            vertical=True
        self.print_aperture(row, vertical)
        self.print_xsec_bias('')
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
        self.print_xsec_bias('')
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
            else:
                self.file.write(', bias="allBias"')

    def print_field(self, row, ndim):
         self.file.write(row.element+'field: field, type="bmap'+str(ndim)+'d", bScaling=1.0, magneticFile="bdsim'+str(ndim)+'d:../magnet_responses/'+row.element+'.dat"')
         self.endl()
       
    def print_fieldmapgeom(self, row, ndim):
        self.line.append(row.element)
        self.print_field(row, ndim)
        self.file.write(row.element+': element, geometryFile="gdml:../'+row.element+'.gdml", fieldAll="'+row.element+'field", l='+str(row.length)+'*mm')
        self.print_xsec_bias('')
        self.endl()

    def print_fieldmap(self, row, magtype, ndim):
        self.line.append(row.element)
        self.print_field(row, ndim)
        self.file.write(row.element+': '+magtype+', fieldVacuum="'+row.element+'field", l='+str(row.length)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt))
        if(not geometry or row.geometry == 'none'):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('')
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

            if(row.element == self.terminal_element): #if we only want the first section of the beam 
                break

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
        elif(use_previous_best_fit):
            self.print_beam()
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
            self.file.write('vacBias: xsecBias, particle="proton", proc="all", xsecfact=100, flag=2;\n')
            self.file.write('matBias: xsecBias, particle="proton", proc="all", xsecfact=100, flag=2;\n')
            self.file.write('allBias: xsecBias, particle="proton", proc="all", xsecfact=100, flag=2;\n')

        self.file.close()

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


    magset["arc_dipole"] = 0.1  ##dummy values for a fake magnet

    magset["QFQ1"] = 0.1 
    magset["BFV1"] = 0.1 
    magset["BFH1"] = 0.1 
    magset["BFV2"] = 0.1 
    magset["QFQ2"] = 0.1 
    magset["QFQ3"] = 0.1 
    magset["BFH2"] = 0.1 
    magset["BFVD1"] = 0.1 
    magset["QFQ4"] = 0.1 
    magset["BFVD2"] = 0.1 
    kvals = {}
    
    if(use_previous_best_fit):
        #use a previous fit as the parameters
        file = ROOT.TFile("../bdsim_optimiser/fit_results.root", "READ")
        tree = file.Get("parameters")
        name = ROOT.TString()
        physical_value = array("d", [0])

        tree.SetBranchAddress("name", name)
        tree.SetBranchAddress("physical_value", physical_value)

        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            kvals[name.Data()] = physical_value[0]

        file.Close()

    else:
        for magnet in magset:
            mag_df = magnet_response[magnet_response['element'] == magnet]
            kvals[magnet] = np.interp(magset[magnet], mag_df['current'], mag_df['kval'])
            zero_field = np.interp(0, mag_df['current'], mag_df['kval'])
            if magnet[0] == 'B': #bending magnets
                kvals[magnet] = -(kvals[magnet]-zero_field) * (proton_momentum/0.2998)  / (0.001*nom.line.loc[nom.line['element'] == magnet].iloc[0]['polelength'])
            else:

                print(f"{magnet} with field {(kvals[magnet]-zero_field)}")
                kvals[magnet] = (kvals[magnet]-zero_field) / (0.001*nom.line.loc[nom.line['element'] == magnet].iloc[0]['polelength'])
            if abs(kvals[magnet]) < 1e-3: #if the strength is zero bdsim will treat it as a drift so force it to be non-zero, if its too small the integrator will fall over however
                kvals[magnet] = 1e-3
    print(kvals)
#    nom.draw_beamline_s(1)
#    nom.draw_beamline_s(2)
#    nom.draw_beamline(0, 1)
#    nom.draw_beamline(0, 2)
#    nom.draw_beamline(1, 2)

    #exit(0)
    if(use_previous_best_fit):
        nom.print_beamline(kvals, "optimised.gmad")
    else:
        nom.print_beamline(kvals, "unoptimised.gmad")
   


#    nom.draw_beamline_s(1)
#    nom.draw_beamline_s(2)
#    nom.draw_beamline(0, 1)
#    nom.draw_beamline(0, 2)

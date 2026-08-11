from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import copy as cp
import sys
import ast
import ROOT
from array import array 
#sys.path.append('../')
#matplotlib.use('qtagg')
pd.set_option('display.max_rows', 500000)
#pd.set_option('display.max_columns', 500000)
false=False
true=True

#if you want to see a plot of the positions of the magnets compared to the survey
plot_beamline=False

#the last element in the beamline you want to include, will place a dump immediately after this
#for the whole beamline set to t2k_target
terminal_element = 'bellows_super'

#known constants
proton_momentum = 30.924 # momentum for a 30GeV KE proton 
vacuum_pressure = 1e-9 #vacuum pressure in bar

ssem_in=True

fit_configuration = False
beam_orbit = False
beam_loss = True
primaries = False

#fit configuration
if(fit_configuration):
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

#beam orbit plot configuration
if(beam_orbit):
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

#beam loss configuration
if(beam_loss):
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
if(primaries):
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
    if(s == 'NA'):
        return []
    s = s.replace(' ', ',')
    s = s.replace('][', '],[')
    elements = ast.literal_eval(s)
    return elements

def string_to_string_list(s):
    if(s == 'NA'):
        return None
    s = s.replace(' ', ',')
    s = s.replace('][', '],[')
    retval = []
    for element in [item.strip('[]') for item in s.split('],[')]:
        elements = element.split(',')
        quoted_elements = [e for e in elements]
        retval.append(quoted_elements)

#    print("string to string list ", retval)
    return retval

def string_to_bool(s):
    if(s == 'True' or s == 'true'):
        return True
    if(s == 'False' or s == 'false'):
        return False
    raise ValueError(f"Invalid string for boolean conversion: {s}")


def read_excel(filename):
    surv = pd.ExcelFile(filename)
    surv = surv.parse('基準点成果表 (側壁fit) ', skiprows=3)
    surv.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters
    surv.columns = ['index', 'ID', 'name', 'x2022', 'y2022', 'h2022', 'x2017', 'y2017', 'h2017', 'x2014', 'y2014', 'h2014', 'dx', 'dy', 'dh', 'comment']
    surv = surv.drop(['x2014', 'y2014', 'h2014'], axis=1)
    surv = surv.drop(['x2017', 'y2017', 'h2017'], axis=1)
    surv = surv.drop(['dx', 'dy', 'dh'], axis=1)
    return surv


class DrawingBeamline:

    def __init__(self):

        self.neutrino_elements  = {"BPV1": [600445.8542, 661508.4667, 0.0], #from Neutrino.20210608.nishida.dwg using the centers of the blue crosshairs
                                   "BPH2": [595809.7446, 661149.3152, 0.0],
                                   "QPQ1": [592868.5575, 660921.4667, 0.0],
                                   "QPQ2": [588581.4019, 660589.3481, 0.0],
                                   "BPD1": [584446.2440, 660250.1661, 0.0],
                                   "BPD2": [580272.1553, 659786.2377, 0.0],
                                   "QPQ3": [569778.2487, 658282.8420, 0.0],
                                   "BPV2": [565423.8316, 657651.1374, 0.0],
                                   "QPQ4": [562365.8431, 657207.5084, 0.0],
                                   "BPH3": [554884.1627, 656122.1248, 0.0],
                                   "QPQ5": [552607.9900, 655791.9158, 0.0],
                                   "BAD1": [545999.7266, 654802.0648, 0.0],
                                   "BAD6": [498651.3247, 633423.4405, 0.0],
                                   "QFQ1": [459488.4326, 551102.2680, 0.0],
                                   "BFV1": [459441.3292, 548401.9382, 0.0],
                                   "BFH1": [459404.7534, 546302.3325, 0.0],
                                   "BFV2": [459300.8919, 540343.1619, 0.0],
                                   "QFQ2": [459249.4660, 537393.6121, 0.0],
                                   "QFQ3": [459178.1699, 532644.1315, 0.0],
                                   "BFH2": [459113.5602, 529594.7944, 0.0],
                                   "BFVD1": [459071.3901, 527175.2420, 0.0],
                                   "QFQ4": [459023.4792, 524426.2999, 0.0]}


        self.prep_elements = {"BPV1": [24863.57600486618,10359.963664168285,0.0], #info from PV1-PQ5 dwg file
                         "BPH2": [20217.26170091104,10545.028987369675,0.0],
                         "QPQ1": [17269.596796088466,10662.380978673755,0.0],
                         "QPQ2": [12973.00241423364,10833.485387679655,0.0],
                         "BPD1": [8825.428297907107,10998.633519809764,0.0],
                         "BPD2": [4626.034704560125,11025.057932632817,0.0],
                         "QPQ3": [-5970.5436311655685,10736.572677916973,0.0],
                         "BPV2": [-10368.885537586993,10616.827233954982,0.0],
                         "QPQ4": [-13457.769351769137,10532.731960884865,0.0],
                         "BPH3": [-21014.86764622477,10327.116821280746,0.0],
                         "QPQ5": [-23313.99341249121,10264.394948653753,0.0], 
                         "BAD1": [-30002.7611,10073.1568,0.0], #first SC arc element defocus 
                         "BAF1": [-34390.468, 9737.2323, 0.0] #second SC arc element focus
                         }
        self.arc_elements = {"BPH2": [595809.8067, 661143.2853, 0.0], #need some way to match the different drawings coordinates
                             "BAD1": [545991.9260, 654791.6516, 0.0], #from Neutrino.061228.ichikawa.dwg
                             "BAF1": [541674.5401, 653945.8497, 0.0],
                             "BAD2": [535756.7003, 652468.4052, 0.0],
                             "BAF2": [531545.9861, 651193.5439, 0.0],
                             "BAD3": [525806.3081, 649129.6952, 0.0],
                             "BAF3": [521744.8388, 647438.6215, 0.0],
                             "BAD4": [516241.2450, 644809.1754, 0.0],
                             "BAF4": [512370.0141, 642719.0779, 0.0],
                             "BAD5": [507158.1433, 639550.4710, 0.0],
                             "BAF5": [503516.1699, 637082.4096, 0.0],
                             "BAD6": [498648.6388, 633406.7793, 0.0],
                             "BAF6": [495272.7198, 630585.6828, 0.0],
                             "BAD7": [490798.6643, 626440.0579, 0.0],
                             "BAF7": [487722.9681, 623294.3735, 0.0],
                             "BAD8": [483687.5221, 618720.6787, 0.0],
                             "BAF8": [480943.0958, 615282.1606, 0.0],
                             "BAD9": [477387.0520, 610326.5739, 0.0],
                             "BAF9": [475001.5373, 606630.0061, 0.0],
                             "BAD10": [471960.8117, 601342.5198, 0.0],
                             "BAF10": [469958.3506, 597425.2065, 0.0],
                             "BAD11": [467463.6339, 591859.2572, 0.0],
                             "BAF11": [465864.3999, 587760.7565, 0.0],
                             "BAD12": [463940.8961, 581972.5336, 0.0],
                             "BAF12": [462761.0979, 577734.1977, 0.0],
                             "BAD13": [461428.2090, 571782.1620, 0.0],
                             "BAF13": [460679.7489, 567446.8454, 0.0],
                             "BAD14": [459950.9203, 561391.1116, 0.0],
                             "BAF14": [459641.3327, 557002.5349, 0.0],
                             "BFH1": [459405.5527, 546306.2492, 0.0], #again, for aligning the different drawings
                             "QFQ3": [(459186.9995+459143.0457)/2.0, (533898.2211+531398.5311)/2.0, 0.0]}
        
        self.ff_elements = {"QFQ1": [24145.1391, 41592.1021, 0.0], #from TSDVBD.20210527.nishida.dwg
                            "BFV1": [21445.1396, 41588.2363, 0.0],
                            "BFH1": [19345.4030, 41585.2070, 0.0],
                            "BFV2": [13385.1431, 41585.2001, 0.0],
                            "QFQ2": [10435.1563, 41585.2104, 0.0],
                            "QFQ3": [5685.0487, 41585.2676, 0.0],
                            "BFH2": [2635.1506, 41585.2143, 0.0],
                            "BFVD1": [215.2366, 41585.2088, -6.31],
                            "QFQ4": [-2534.1466, 41585.2222, -69.60],
                            "BFVD2": [-6082.1474, 41585.2100, -173.7408], ##OLD FVD2 position!!
                            "NUTGT": [-12673.4740, 41585.2097, -570.2937]}

        #convert these elements to numpy arrays for easier manipulation
        for element in self.neutrino_elements:
            self.neutrino_elements[element] = np.array(self.neutrino_elements[element])
        for element in self.prep_elements:
            self.prep_elements[element] = np.array(self.prep_elements[element])
        for element in self.arc_elements:
            self.arc_elements[element] = np.array(self.arc_elements[element])
        for element in self.ff_elements:
            self.ff_elements[element] = np.array(self.ff_elements[element])

        self.align_elements()
        self.drop_redundant_elements()


    def align_elements(self):
        #first let's set use QPQ1 as the origin
        reference_point = cp.deepcopy(self.prep_elements["QPQ1"])
        for key, element in self.prep_elements.items():
            element -= reference_point
        #now rotate so that the beam is along the x axis (looks more like the drawings that way)
        xvec = self.prep_elements["QPQ2"] - self.prep_elements["QPQ1"]
        angle = np.arctan2(xvec[1], xvec[0])
        R = rotmatxy3d(angle)
        for key, element in self.prep_elements.items():
            #print(key, self.prep_elements[key], "before rotation")
            self.prep_elements[key] = R.dot(self.prep_elements[key])
            #print(key, self.prep_elements[key], "after rotation")
        #add offset of 10450mm so that 0 is at NUIN (this agrees with SAD perfectly and the drawings to 0.8mm)
        offset = np.array([10450.0, 0.0, 0.0])
        for element in self.prep_elements:
            self.prep_elements[element] += offset

        #now align the vector connecting BPH2 and BAD1 to be the same between the prep and arc elements
        prep_vec = self.prep_elements["BAD1"] - self.prep_elements["BPH2"]
        arc_vec = self.arc_elements["BAD1"] - self.arc_elements["BPH2"]
        angle2 = np.arctan2(arc_vec[1], arc_vec[0]) - np.arctan2(prep_vec[1], prep_vec[0])
        R2 = rotmatxy3d(angle2)
        for element in self.arc_elements:
            self.arc_elements[element] = R2.dot(self.arc_elements[element])

        #now align the position of BPH2 between the two sets
        shift = self.prep_elements["BPH2"] - self.arc_elements["BPH2"]
        for element in self.arc_elements:
            self.arc_elements[element] += shift

        #now align the vector connecting BFH1 and QFQ3 to be the same between the arc and ff elements
        arc_vec2 = self.arc_elements["QFQ3"] - self.arc_elements["BFH1"]
        ff_vec = self.ff_elements["QFQ3"] - self.ff_elements["BFH1"]
        angle3 = np.arctan2(ff_vec[1], ff_vec[0]) - np.arctan2(arc_vec2[1], arc_vec2[0])
        R3 = rotmatxy3d(angle3)
        for element in self.ff_elements:
            self.ff_elements[element] = R3.dot(self.ff_elements[element])

        #now align the position of BFH1 between the two sets
        shift2 = self.arc_elements["BFH1"] - self.ff_elements["BFH1"]
        for element in self.ff_elements:
            self.ff_elements[element] += shift2
        #should now have a consistent coordinate system across the sections with the entry to the neutrino beamline
        #being along the [1, 0, 0] axis and all elements being at height=0


        #now for the other version of the drawings
        vec = self.neutrino_elements["QPQ2"] - self.neutrino_elements["QPQ1"]
        angle = np.arctan2(vec[1], vec[0])
        R4 = rotmatxy3d(angle)
        for element in self.neutrino_elements:
            self.neutrino_elements[element] = R4.dot(self.neutrino_elements[element])
        #align QPQ1
        shift = self.neutrino_elements["QPQ1"] - np.array([10450.0, 0.0, 0.0])
        for element in self.neutrino_elements:
            self.neutrino_elements[element] -= shift

        print(f"alternative drawing beamline positions: {self.neutrino_elements}")


    def drop_redundant_elements(self):
        #in theory we now have a single coherent coordinate system for all elements
        #so lets drop the repeated elements
        del self.prep_elements["BAD1"]
        del self.prep_elements["BAF1"]
        del self.arc_elements["BPH2"]
        del self.arc_elements["QFQ3"]
        del self.arc_elements["BFH1"]

        #and combine these dicts into a single beamline dict
        #this is the 3d position of the magnet centers according to the various beamline design drawings
        self.beamline = {**self.prep_elements, **self.arc_elements, **self.ff_elements}
        #remove the different sections to avoid confusion, from now on, beamline is the only thing that matters
        del self.prep_elements
        del self.arc_elements
        del self.ff_elements
        print("Final aligned beamline positions from drawings: ", self.beamline)


    def plot_drawing_beamline(self, x, y, plot=None):
        labels = ["X (mm)", "Y (mm)", "Z (mm)"]
        if(plot is None):
            plt.figure()
        plt.plot(-999., -999., 'rs', alpha=0.5, label='Beamline Drawings (wrong)')
        plt.plot(-999., -999., 'bs', alpha=0.5, label='Beamline Drawings (more correct? but incomplete)')
        if(self.beamline):
            for element in self.beamline.keys():
                plt.plot(self.beamline[element][x], self.beamline[element][y], 'rs', alpha=0.5)
                plt.text(self.beamline[element][x], self.beamline[element][y], element, rotation=-45, horizontalalignment='right')
        else:
            for element in self.prep_elements:
                plt.plot(self.prep_elements[element][x], self.prep_elements[element][y], 'rs', alpha=0.5)
                plt.text(self.prep_elements[element][x], self.prep_elements[element][y], element, rotation=-45, horizontalalignment='right')
            for element in self.arc_elements:
                plt.plot(self.arc_elements[element][x], self.arc_elements[element][y], 'rs', alpha=0.5)
                plt.text(self.arc_elements[element][x], self.arc_elements[element][y], element, rotation=-45, horizontalalignment='right')
            for element in self.ff_elements:
                plt.plot(self.ff_elements[element][x], self.ff_elements[element][y], 'rs', alpha=0.5)
                plt.text(self.ff_elements[element][x], self.ff_elements[element][y], element, rotation=-45, horizontalalignment='right')
        for element in self.neutrino_elements:
            plt.plot(self.neutrino_elements[element][x], self.neutrino_elements[element][y], 'bs', alpha=0.5)
            plt.text(self.neutrino_elements[element][x], self.neutrino_elements[element][y], element, rotation=-45, horizontalalignment='right')
        
        plt.xlabel(labels[x])
        plt.ylabel(labels[y])
        plt.title('Aligned Beamline Elements')
#        plt.axis('equal')
        plt.grid()
        if(plot is None):
            plt.legend()
            plt.show()

class SpreadsheetBeamline:

    def __init__(self, file):
        self.line = strip_whitespace(pd.read_csv(file, header=0, skipinitialspace=True)) #the csv containing beampipe properties
        self.line['s_start'] = self.line['length'].shift().cumsum()
        self.line['center_pos'] = None
        self.line['misalign'] = None

        #apply the function to convert the survey offset strings to lists
        self.line['survey_offset'] = self.line['survey_offset'].apply(string_to_list)
        self.line['offset'] = self.line['offset'].apply(string_to_list)
        self.line['survey_name'] = self.line['survey_name'].apply(string_to_string_list)
        self.line['geometry'] = self.line['geometry'].apply(string_to_bool)
        self.line['mark'] = pd.to_numeric(self.line['mark'], errors='coerce')
        self.line['polelength'] = pd.to_numeric(self.line['polelength'], errors='coerce')
        self.line['beamline_dir'] = None
        self.line = self.line.set_index('element', drop=False)

        #if the element 'Survey_offset' is not NaN then this element has a valid survey measurement
        self.survey_elements = self.line[self.line['survey_name'].notna()]

        self.calculate_spreadsheet_beamline()

    def get_bellows(self, l, offset, angle):
        alpha = 2.0*np.arcsin(offset/l)
        angle = angle / 2.0
        bend1 = angle - alpha
        bend2 = angle + alpha
        return bend1, bend2


    def calculate_spreadsheet_beamline(self):
        current_s = 0.
        current_pos = np.array([0., 0., 0.])
        current_dir = np.array([1., 0., 0.])

        #we will keep track of the acumulated bend angles
        #NOTE this only works if there are no horisontal bends after vertical bends
        #this is true for the T2K beamline but not true in general
        #we'll keep a vector based on first applying the acumulated horisontal bends and then
        #applying the acumulated vertical bends in current_dir
        acum_horisontal_bend = 0.
        acum_vertical_bend = 0.
        for index, element in self.line.iterrows():
            if element['type'] == 'drift' or element['type'] == 'ssem' or element['type'] == 'wsem':
                current_pos += current_dir * element['length']
                current_s += element['length']
                self.line.at[index, 'beamline_dir'] = cp.deepcopy(current_dir)

            elif element['type'] == 'bellows':
                #these allow us to align bending magnets with the survey data, essentially make the combination of a 
                #'c' and 's' shaped duct, 'c' for half the bending angle of the magnet and 's' for the misalignment of 
                #the magnet center in to the survey. It is composed of two rbend elements that together encompass
                #the bend angle desired and the positional shift relative to the previous bit of beamline 
                #it really is like the bellows, allowing for flexibility in angle and position between elements

                #first calculate the 'c' part, half of half of the bend angle of the magnet since we're using two rbends for the bellows
                #every magnet has a bellows before and after (pretty much) so assign half to each bellows (this is done in the spreadsheet)

                bend1, bend2 = self.get_bellows(element['length'], element['offset'][0], element['angle'])
                if(element['tilt'] > 1): #assume this is a vertical tilting bellows (should only be BFVD1 and BFVD2)
                    acum_vertical_bend += bend1/2.0
                else: #horisontal bend
                    acum_horisontal_bend += bend1/2.0

                current_dir = np.array([np.cos(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_vertical_bend)])
                current_pos += current_dir * element['length']/4.0
                current_pos += current_dir * element['length']/4.0

                if(element['tilt'] > 1): #assume this is a vertical tilting bellows (should only be BFVD1 and BFVD2)
                    acum_vertical_bend += bend1/2.0
                else: #horisontal bend
                    acum_horisontal_bend += bend1/2.0


                current_dir = np.array([np.cos(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_vertical_bend)])

                if(element['tilt'] > 1): #assume this is a vertical tilting bellows (should only be BFVD1 and BFVD2)
                    acum_vertical_bend += bend2/2.0
                else: #horisontal bend
                    acum_horisontal_bend += bend2/2.0


                current_dir = np.array([np.cos(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_vertical_bend)])
                current_pos += current_dir * element['length']/4.0
                current_pos += current_dir * element['length']/4.0

                if(element['tilt'] > 1): #assume this is a vertical tilting bellows (should only be BFVD1 and BFVD2)
                    acum_vertical_bend += bend2/2.0
                else: #horisontal bend
                    acum_horisontal_bend += bend2/2.0

                current_dir = np.array([np.cos(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_vertical_bend)])
                self.line.at[index, 'beamline_dir'] = cp.deepcopy(current_dir)
                current_s += element['length']

            elif element['type'] == 'rbend' or element['type'] == 'sbend':
                #first apply the upstream drift
                drift_u = element['mark'] - element['polelength']/2.
                current_pos += current_dir * drift_u
                current_s += drift_u

                bend_angle = element['angle']
                bend_fudge_factor = -0.0e-5 #3.097e-4/2.#1.42e-4 #radians, to make the spreadsheet and survey line up better
                if(element['element'].startswith('BAF')):
                   bend_angle = bend_angle + bend_fudge_factor #apply fudge factor to arc magnets only
                if(element['element'].startswith('BAD')):
                   bend_angle = bend_angle - bend_fudge_factor #apply fudge factor to arc magnets only
                #now apply half of the bend angle
                if element['tilt'] > 1: #assume this is a vertical tilting magnet
                    acum_vertical_bend += bend_angle/2.
                else: #horisontal bend
                    acum_horisontal_bend += bend_angle/2.
                
                current_dir = np.array([np.cos(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_vertical_bend)])
                current_pos += current_dir * element['polelength']/2.

                self.line.at[index, 'beamline_dir'] = cp.deepcopy(current_dir)

                current_s += element['polelength']/2.
                #save magnet center position
                self.line.at[index, 'center_pos'] = cp.deepcopy(current_pos)

                #rest of the magnet
                current_pos += current_dir * element['polelength']/2.
                current_s += element['polelength']/2.


                #now apply the downstream drift
                #again half of the bend angle
                if element['tilt'] > 1: #assume this is a vertical tilting magnet
                    acum_vertical_bend += bend_angle/2.
                else: #horisontal bend
                    acum_horisontal_bend += bend_angle/2.
                
                current_dir = np.array([np.cos(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_horisontal_bend)*np.cos(acum_vertical_bend), -np.sin(acum_vertical_bend)])
                drift_d = element['length'] - element['mark'] - element['polelength']/2.

                drift_fudge_factor = 0.0
                if(element['element'].startswith('BA')):
                   drift_d += drift_fudge_factor #apply fudge factor to arc magnets only
                current_pos += current_dir * drift_d
                current_s += drift_d

            elif element['type'] == 'quadrupole':
                #essentially two drifts in a row
                drift_u = element['mark']
                current_pos += current_dir * drift_u
                current_s += drift_u

                self.line.at[index, 'beamline_dir'] = cp.deepcopy(current_dir)
                #save magnet center position
                self.line.at[index, 'center_pos'] = cp.deepcopy(current_pos)


                #now the second drift
                drift_d = element['length'] - element['mark']
                current_pos += current_dir * drift_d
                current_s += drift_d
            elif element['type'] == 'dump':
                self.line.at[index, 'center_pos'] = cp.deepcopy(current_pos)
                self.line.at[index, 'beamline_dir'] = cp.deepcopy(current_dir)
            else:
                print(f"WARNING: unhandled element type {element['type']} in nominal beamline calculation")

    def translate_spreadsheet_beamline(self, shift):
        for indx, element in self.line.iterrows():
            if element['center_pos'] is not None:
                element['center_pos'] += shift

    def rotate_spreadsheet_beamline(self, position, R):
        self.translate_spreadsheet_beamline(-position) #shift origin to 'position'
        for idx, element in self.line.iterrows():
            if element['center_pos'] is not None:
                element['center_pos'] = R.dot(element['center_pos'])
        self.translate_spreadsheet_beamline(position) #shift back to original origin

    def plot_spreadsheet_beamline(self, x, y, plot=None):
        labels = ["X (mm)", "Y (mm)", "Z (mm)"]
        if(plot is None):
            plt.figure()
        plt.plot(-999., -999., 'go', alpha=0.5, label='Fujii-Sans spreadsheet')
        for index, element in self.line.iterrows():
            if(element['center_pos'] is None):
                continue
            plt.plot(element['center_pos'][x], element['center_pos'][y], 'go', alpha=0.5)
            plt.text(element['center_pos'][x], element['center_pos'][y], element['element'], rotation=-45, horizontalalignment='right')
            plt.arrow(element['center_pos'][x], element['center_pos'][y], 100*element['beamline_dir'][x], 100*element['beamline_dir'][y], head_width=2, head_length=10, fc='g', ec='g', alpha=0.5)

        plt.xlabel(labels[x])
        plt.ylabel(labels[y])
        plt.title('Beamline Elements from Spreadsheet')
#        plt.axis('equal')
        plt.grid()
        if(plot is None):
            plt.show()


class SurveyBeamline:
    survey_data = {}
    survey_centers = {}

    def __init__(self, drawing_beamline, spreadsheet_beamline):
        self.survey = read_excel("03_2022_Neutrino.xlsx")
        self.drawing_beamline = drawing_beamline
        self.spreadsheet_beamline = spreadsheet_beamline
        self.elements = cp.deepcopy(spreadsheet_beamline.survey_elements)
        self.load_survey()

    def load_survey(self):
        for index, element in self.elements.iterrows():
#            print(element['element'])
            self.survey_data[element['element']] = []
            for point in element['survey_name'][0]:
#                print("survey point ", point)
                self.survey_data[element['element']].append(np.array([self.survey[self.survey['name'] == point]['x2022'].values,
                                            self.survey[self.survey['name'] == point]['y2022'].values,
                                            self.survey[self.survey['name'] == point]['h2022'].values]))
                                                                      
            self.survey_data[element['element']] = np.array(self.survey_data[element['element']])
            self.survey_data[element['element']] = self.survey_data[element['element']].squeeze()

        print("survey data read in:", self.survey_data)
        for index, row in self.survey.iterrows():
            name = row['name']
            x = row['x2022']
            y = row['y2022']
            h = row['h2022']

    def initial_alignment(self):
        #a rough alignment of the survey points with the nominal beamline (can only be rough since we next need to calculate magnet centers)

        PV1 = self.survey_data['BPV1']
        PV1_vertical_offset = np.array(self.spreadsheet_beamline.survey_elements[self.spreadsheet_beamline.survey_elements['element'] == 'BPV1']['survey_offset'].tolist()).squeeze()
        PV1_center = PV1.mean(axis=0) -  PV1_vertical_offset[0][2] * np.array([0.0, 0.0, 1.0])

        PH2 = self.survey_data['BPH2']
        PH2_vertical_offset = np.array(self.spreadsheet_beamline.survey_elements[self.spreadsheet_beamline.survey_elements['element'] == 'BPH2']['survey_offset'].tolist()).squeeze()
        PH2_center = PH2.mean(axis=0) -  PH2_vertical_offset[0][2] * np.array([0.0, 0.0, 1.0])

        #first rotate the whole survey so that the line between centers is horisontal and along the x axis
        survey_vec = PH2_center - PV1_center
        thetaxz = np.arctan(survey_vec[2]/np.sqrt(survey_vec[0]**2 + survey_vec[1]**2))
        self.rotate_survey_xz(-thetaxz, self.survey_data)

        thetaxy = np.arctan2(survey_vec[1], survey_vec[0])
        self.rotate_survey_xy(-thetaxy, self.survey_data)
        
        #find expected position of BPV1
        expected_PV1 = self.drawing_beamline.beamline['BPV1']
        PV1_center = self.survey_data['BPV1'].mean(axis=0) -  PV1_vertical_offset[0][2] * np.array([0.0, 0.0, 1.0])
        shift = expected_PV1 - PV1_center
        self.shift_survey(shift, self.survey_data)

        print("After initial alignment:", self.survey_data)

    def precise_alignment(self):
        #will align QPQ1 and QPQ2 to lie in the xz plane
        #then align QPQ1 and QPQ2 to lie in the xy plane
        #the vector QPQ2-QPQ1 is now aligned along [1, 0, 0]
        #extra DOF rotation around z axis
        #rotate yz plane so QFQ3 lies in the xy plane
        #place QPQ1 at its nominal position
        #will then rotate along the [1 0 0] vector to have QFQ3 center in the xy plane
        self.calculate_centers()
        QPQ1 = self.survey_centers['QPQ1']
        QPQ2 = self.survey_centers['QPQ2']
        thetaxy = np.arctan2((QPQ2-QPQ1)[1], (QPQ2-QPQ1)[0])
        print("Precise angle to rotate:", thetaxy)
        self.rotate_survey_xy(-thetaxy, self.survey_data)
        self.calculate_centers()



        QPQ1 = self.survey_centers['QPQ1']
        QPQ5 = self.survey_centers['QPQ5']
        thetaxz = np.arctan2((QPQ5-QPQ1)[2], (QPQ5-QPQ1)[0])

        print("Precise angle to rotate:", thetaxz)
        self.rotate_survey_xz(-thetaxz, self.survey_data)
        self.calculate_centers()

       


        #can now be a bit more precise, QFQ3 is much further from QPQ1 than QPQ2 so give up on QPQ2 center having no horisontal shift and no vertical shift
#        QPQ1 = self.survey_centers['QPQ1']
#        QFQ3 = self.survey_centers['QFQ3']
#
#        spreadsheet_vector = normvec((self.spreadsheet_beamline.line.loc['QFQ3', 'center_pos'] - self.spreadsheet_beamline.line.loc['QFQ3', 'center_pos'])[0:2])
#        design_vector = normvec((self.drawing_beamline.beamline['QFQ3'] - self.drawing_beamline.beamline['QPQ1'])[0:2]) #only care about xy
#        actual_vector = normvec((QFQ3-QPQ1)[0:2])
#        print(f"Aligning QFQ3 with QPQ1, current actual direction: {actual_vector} and from the drawings this should be {design_vector}")
#        angle_diff = np.arccos(actual_vector.dot(design_vector))
#        print(f"angle between QFQ3-QPQ1 actual vs design = {angle_diff} rad")
#        spreadheet_design_angle_diff = np.arccos(spreadsheet_vector.dot(design_vector))
#        R = rotmatxy3d(-spreadheet_design_angle_diff)
#        self.spreadsheet_beamline.rotate_spreadsheet_beamline(cp.deepcopy(self.spreadsheet_beamline.line.loc['QFQ1', 'center_pos']), R)
#        self.rotate_survey_xy(-angle_diff, self.survey_data)
#        self.calculate_centers()


        QPQ1 = self.survey_centers['QPQ1']
        QFQ3 = self.survey_centers['QFQ3']

        thetayz = np.arctan2((QFQ3-QPQ1)[2], (QFQ3-QPQ1)[1])
        print("Precise angle to rotate yz:", thetayz)
        self.rotate_survey_yz(-thetayz, self.survey_data)
        self.calculate_centers()


        QPQ1 = self.survey_centers['QPQ1']
        QPQ2 = self.survey_centers['QPQ2']
        thetaxy = np.arctan2((QPQ2-QPQ1)[1], (QPQ2-QPQ1)[0])
        print("Precise angle to rotate:", thetaxy)
        self.rotate_survey_xy(-thetaxy, self.survey_data)
        self.calculate_centers()


        #now place QPQ1 at sQPQ1*[1, 0, 0]
        QPQ1 = self.survey_centers['QPQ1']
        survey_shift = self.drawing_beamline.beamline['QPQ1'] - QPQ1
        self.shift_survey(survey_shift, self.survey_data)
        self.calculate_centers()

    def calculate_centers(self):
        #this evaluates the magnet center position based on the survey points and nominal offsets
        #a number of simplifying assumptions are used: 
        #that survey points on the same magnet are the same distance above the center of the beam
        #that the magnet is level on the axis perpendicular to the beamline (i.e. no roll)
        #that two survey points are present with the same lateral offset from the center (to allow us to define its axis)
        #these statements are true (to some high precision) for the T2K Neutrino beamline but may not be true for others

        for key, magnet in self.survey_data.items():
#            print(f"Calculating center for {key}")
            #first we need a guess of the long-axis of the magnet
            if(magnet.shape[0] == 2): #a magnet with two survey points
                s = magnet[1] - magnet[0]
                s /= np.linalg.norm(s)
                lateral = np.array([s[1], -s[0], 0.0]) #perpendicular to the beam direction in the horizontal plane
                lateral /= np.linalg.norm(lateral)
            elif(magnet.shape[0] == 4): #a magnet with four survey points (two on each side)
                s = magnet[2] - magnet[0] #grab two points on the same side
                s /= np.linalg.norm(s)
                lateral = magnet[0] - magnet[1] #vector between the two points on the same end
                lateral /= np.linalg.norm(lateral)
            else:
                raise Exception(f"Magnet {key} has invalid number of survey points {magnet.shape[0]}")

            vertical = -np.cross(s, lateral) #it is possible that the magnet has pitch but we assume no roll so vertical is defined by cross product of s and lateral
            vertical /= np.linalg.norm(vertical)
            point_offsets = np.array(self.spreadsheet_beamline.survey_elements[self.spreadsheet_beamline.survey_elements['element'] == key]['survey_offset'].tolist()).squeeze()
            self.survey_centers[key] = magnet[0] - point_offsets[0][0]*s - point_offsets[0][1]*lateral - point_offsets[0][2]*vertical;

#        self.print_survey_vs_drawing_centers("QPQ1")
#        self.print_survey("BFVD2")
#        self.print_survey_vs_drawing_centers("QFQ4")

    def print_survey(self, key=None):
        if(key is not None):
            print(f"Survey points for {key} are {self.survey_data[key]}")
        else:
            for key in self.survey_data:
                self.print_survey(key)

    def print_survey_vs_drawing_centers(self, key=None):
        if(key is not None):
            print(f"Survey center for {key} is {self.survey_centers[key]}")
            print(f"Nominal center for {key} is {self.drawing_beamline.beamline[key]}")
            print(f"Difference for {key} is {self.survey_centers[key]-self.drawing_beamline.beamline[key]}")
        else:
            for key in self.survey_centers:
                self.print_survey_vs_drawing_centers(key)

    def calculate_spreadsheet_misalignments(self):
            for key in self.survey_centers:
                difference = self.survey_centers[key]-self.spreadsheet_beamline.line.loc[key, 'center_pos']
                displacement_s = normvec(self.spreadsheet_beamline.line.loc[key, 'beamline_dir']).dot(difference)
                displacement_vertical = difference[2]
                displacement_lateral = difference[0:2] - displacement_s * normvec(self.spreadsheet_beamline.line.loc[key, 'beamline_dir'][0:2])
                displacement_lateral = np.linalg.norm(displacement_lateral)
                self.spreadsheet_beamline.line.at[key, 'misalign'] = cp.deepcopy(np.array([displacement_s, displacement_lateral, displacement_vertical]))

    def print_survey_vs_spreadsheet_centers(self, key=None):
        if(key is not None):
            print(f"Survey center for {key} is {self.survey_centers[key]}")
            print(f"Spreadsheet center for {key} is {self.spreadsheet_beamline.line.loc[key, 'center_pos']}")
            print(f"Difference for {key} is {self.survey_centers[key]-self.spreadsheet_beamline.line.loc[key, 'center_pos']}")
            print(f"Misalignment for {key} is [s, horisontal, vertical]: {self.spreadsheet_beamline.line.loc[key, 'misalign']}")
        else:
            for key in self.survey_centers:
                self.print_survey_vs_spreadsheet_centers(key)
                print("")


    def rotate_survey_xz(self, angle, data):
        for element in data:
            points = data[element]
            c = np.cos(angle)
            s = np.sin(angle)
            R = np.array([[c, 0, -s],
                          [0, 1, 0],
                          [s, 0, c]])
            rotated_points = R.dot(points.T).T
            data[element] = rotated_points

        
    def rotate_survey_xy(self, angle, data):
        for element in data:
            points = data[element]
            c = np.cos(angle)
            s = np.sin(angle)
            R = np.array([[c, -s, 0],
                          [s, c, 0],
                          [0, 0, 1]])
            rotated_points = R.dot(points.T).T
            data[element] = rotated_points

    def rotate_survey_yz(self, angle, data):
        for element in data:
            points = data[element]
            c = np.cos(angle)
            s = np.sin(angle)
            R = np.array([[1, 0, 0],
                          [0, c, -s],
                          [0, s, c]])
            rotated_points = R.dot(points.T).T
            data[element] = rotated_points


    def shift_survey(self, shift, data):
        for element in data:
            data[element] = data[element] + shift
    
    def plot_survey(self, x, y, plot=None):
        labels = ["X (mm)", "Y (mm)", "Z (mm)"]
        if(plot is None):
            plt.figure()
        plt.plot(-999., -999., 'bo', alpha=0.5, label='Survey Points')
        for element in self.survey_data:
            plt.plot(self.survey_data[element][:,x], self.survey_data[element][:,y], 'bo', alpha=0.5)
#            plt.text(self.survey_data[element][:,x].mean(), self.survey_data[element][:,y].mean(), element, rotation=-45, horizontalalignment='right')
        try:  #if this has been calculated plot it as well

            plt.plot(-999., -999., 'ro', alpha=0.8, label='Survey Center')
            for center in self.survey_centers:
                plt.plot(self.survey_centers[center][x], self.survey_centers[center][y], 'ro', alpha=0.8)
                plt.text(self.survey_centers[center][x], self.survey_centers[center][y], center, rotation=-45, horizontalalignment='right')
        except:
            pass
        plt.xlabel(labels[x])
        plt.ylabel(labels[y])
        plt.title('Surveyed Beamline Elements')
        #plt.axis('equal')
        plt.grid()
        if(plot is None):
            plt.show()


class BeamlinePrinter:

    file = None
    def __init__(self, beamline, kv, primaries_only=False):
        self.beamline = beamline
        self.kvals = kv
        self.s = 0
        self.blmID = 1
        self.primaries_only=primaries_only
        self.line = []
        self.terminal_element = terminal_element

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

      X0={self.kvals["X0"]}*m,
      Xp0={self.kvals["Xp0"]},
      emitx={self.kvals["emitx"]}*m*rad,
      betx={self.kvals["betx"]}*m,
      alfx={self.kvals["alfx"]},
      dispx=0.423734*m,
      dispxp=0.0719639,

      Y0={self.kvals["Y0"]}*m,
      Yp0={self.kvals["Yp0"]},
      emity={self.kvals["emity"]}*m*rad,
      bety={self.kvals["bety"]}*m,
      alfy={self.kvals["alfy"]},
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

    def print_drift(self, row, name, driftlen, misalign=[0., 0., 0.]):
        self.line.append(name)
        self.file.write(name + ': drift, l=' + str(driftlen)+ '*mm, offsetX='+str(misalign[1])+'*mm, offsetY='+str(misalign[2])+'*mm')
        self.print_aperture(row)
        self.print_xsec_bias('')

    def print_bellows(self, row, angle1, angle2):
        name = row.element + "_ubellows"
        self.line.append(name)
        self.file.write(name + ': rbend, l=' + str(row.length/2.0) + '*mm, angle=' + str(angle1) + ', tilt=' + str(row.tilt) + ', B=0.001*T, magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('')
        self.endl()

        name = row.element + "_dbellows"
        self.line.append(name)
        self.file.write(name + ': rbend, l=' + str(row.length/2.0) + '*mm, angle=' + str(angle2) + ', tilt=' + str(row.tilt) + ', B=-0.001*T, magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('')
        self.endl()

    def print_bend_magnet(self, row):
        self.line.append(row.element)
        if(type(self.kvals[row.element]) == list): #combined function bend + quad
            self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt)+', B='+str(self.kvals[row.element][0])+'*T, k1='+str(self.kvals[row.element][1]))
        else: #standard bend
            self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt)+', B='+str(self.kvals[row.element])+'*T')

        if(misalignments and row.angle == 0.0): #bend magnets with actual deflection are handled by the bellows method since BDSIM would collide geometries
            self.file.write(', offsetX='+str(row.misalign[1])+'*mm, offsetY='+str(row.misalign[2])+'*mm')
        if(not geometry or row.geometry==False):
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
        if(not geometry or row.geometry == False):
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
        #CERN TODO
        misalign = [0., 0., 0.]
#        misalign = row.misalign
#        if(np.isnan(row.misalign).any()):
#            misalign = [0., 0., 0.]
        if(ssem_in):
            self.print_target(row.element, str(thickness)+"*mm", "G4_Ti", row.aperture_x, misalign)
        else:
            self.print_drift(row, row.element, thickness)
            self.endl()
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
        if(not geometry or row.geometry == False):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('')
        self.endl()

    def print_blms(self, row):
        dx = string_to_list(row.blm_offset_x)
        dy = string_to_list(row.blm_offset_y)
        ds = [str(float(entry)+row.polelength/2.0) for entry in string_to_list(row.blm_offset_s)]  #TODO THIS ONLY WORKS FOR PURE QUAD/DIPOLE FIELDMAP WILL BREAK THIS
        orientation = string_to_string_list(row.blm_orientation)[0]

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

    def print(self, filename):
        self.file = open(filename, "w")
        self.file.write('chrg: scorer, type="cellcharge";\n')
        self.file.write('eDep: scorer, type="depositeddose";\n')
        previously_drift = False
        prev_row = []
        driftlen = 0.0
        for row in self.beamline.spreadsheet_beamline.line.itertuples():
            if(row.type != 'drift' and previously_drift and merge_drifts):
                #we have reached a non-drift element, so flush the drift components to the file
                self.print_drift(prev_row, prev_row.element, driftlen)
                self.endl()
                driftlen = 0.0
            previously_drift = False
            if(row.type in ['rbend', 'sbend', 'quadrupole']):
                self.file.write('\n') #give some breathing room
                if(misalignments):
                    self.print_drift(row, row.element+'_driftu', float(row.mark)-row.polelength/2.+row.misalign[0], row.misalign)
                else:
                    self.print_drift(row, row.element+'_driftu', float(row.mark)-row.polelength/2.)

                self.endl()
                if(row.type == 'quadrupole'):
                    self.print_quad_magnet(row)
                else:
                    self.print_bend_magnet(row)
                if(misalignments):
                    self.print_drift(row, row.element+'_driftd', row.length - (float(row.mark)+row.polelength/2.)-row.misalign[0], row.misalign)
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

            elif(row.type == 'bellows'):
                bend1, bend2 = self.beamline.spreadsheet_beamline.get_bellows(row.length, row.offset[0], row.angle)
                self.print_bellows(row, bend1, bend2)
            elif(row.type == 'ssem'):
                self.print_ssem(row, 15e-3)
            elif(row.type == 'wsem'):
                self.print_ssem(row, 1e-4)
            elif(row.type == 'dump'):
                self.print_dump(row)
            else: #otherwise, i.e. drift or collimator
                previously_drift = True
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
                if(row.type != "dump"):
                    dump = row._replace(element="d1", length=1000)
                    self.print_dump(dump)
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
        self.file.write('option, integratorSet="bdsimtwo";\n')
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


    drawing = DrawingBeamline()
    spreadsheet = SpreadsheetBeamline("../fujii-san_hand_tuned.csv")
#just pass the keys for the magnets we want to include in the survey
    surv = SurveyBeamline(drawing, spreadsheet)
    surv.initial_alignment()
    surv.precise_alignment()
    surv.calculate_spreadsheet_misalignments()

    surv.print_survey_vs_spreadsheet_centers()

    fujiisan_ps_dir = normvec(spreadsheet.line.at['QPQ2', 'center_pos'] - spreadsheet.line.loc['QPQ1', 'center_pos'])
    drawing_ps_dir = normvec(drawing.beamline['QPQ2'] - drawing.beamline['QPQ1'])
    survey_ps_dir = normvec(surv.survey_centers['QPQ2'] - surv.survey_centers['QPQ1'])

    fujiisan_ps_after_bpd2 = normvec(spreadsheet.line.at['QPQ5', 'center_pos'] - spreadsheet.line.loc['QPQ3', 'center_pos'])
    survey_ps_after_bpd2 = normvec(surv.survey_centers['QPQ5'] - surv.survey_centers['QPQ3'])

    fujiisan_after_arc = normvec(spreadsheet.line.at['BFV1', 'center_pos'] - spreadsheet.line.loc['QFQ1', 'center_pos'])
    survey_after_arc = normvec(surv.survey_centers['BFV1'] - surv.survey_centers['QFQ1'])


    fujiisan_ff_dir = normvec(spreadsheet.line.at['QFQ3', 'center_pos'] - spreadsheet.line.loc['QFQ2', 'center_pos'])
    drawing_ff_dir = normvec(drawing.beamline['QFQ3'] - drawing.beamline['QFQ2'])
    survey_ff_dir = normvec(surv.survey_centers['QFQ3'] - surv.survey_centers['QFQ2'])

    print(f"FF According to Fujii-san spreadsheet {fujiisan_ff_dir}")
    print(f"FF According to beamline drawings     {drawing_ff_dir}")
    print(f"FF According to survey                {survey_ff_dir}")
    print(f"FF According to Fujii-san spreadsheet {np.arctan(fujiisan_ff_dir[1] / fujiisan_ff_dir[0])}")
    print(f"FF According to beamline drawings     {np.arctan(drawing_ff_dir[1] / drawing_ff_dir[0])}")
    print(f"FF According to survey                {np.arctan(survey_ff_dir[1] / survey_ff_dir[0])}")

    print(f"Bend from PS to after BPD2 Fujii-san spreadsheet {np.arccos(fujiisan_ps_after_bpd2.dot(fujiisan_ps_dir))}")
    print(f"Bend from PS to after BPD2 survey                {np.arccos(survey_ps_after_bpd2.dot(survey_ps_dir))}")

    print(f"Bend from PS to after arc Fujii-san spreadsheet {np.arccos(fujiisan_after_arc.dot(fujiisan_ps_dir))}")
    print(f"Bend from PS to after arc survey                {np.arccos(survey_after_arc.dot(survey_ps_dir))}")


    print(f"Total x bend Fujii-san spreadsheet {np.arccos(fujiisan_ff_dir.dot(fujiisan_ps_dir))}")
    print(f"Total x bend beamline drawings     {np.arccos(drawing_ff_dir.dot(drawing_ps_dir))}")
    print(f"Total x bend survey                {np.arccos(survey_ff_dir.dot(survey_ps_dir))}")


    if(plot_beamline):
        plt.plot()
        drawing.plot_drawing_beamline(0, 1, True)
        spreadsheet.plot_spreadsheet_beamline(0, 1, True)
        surv.plot_survey(0, 1, True)
    #    plt.axis('equal')
        plt.legend()
        plt.show()
        plt.plot()
        drawing.plot_drawing_beamline(1, 2, True)
        spreadsheet.plot_spreadsheet_beamline(1, 2, True)
        surv.plot_survey(1, 2, True)
        plt.legend()
        plt.show()
        plt.plot()
        drawing.plot_drawing_beamline(0, 2, True)
        spreadsheet.plot_spreadsheet_beamline(0, 2, True)
        surv.plot_survey(0, 2, True)
        plt.legend()
        plt.show()



    #run 910216
#    vec_magset = [0 ,
#    -15 ,
#    520 ,
#    0 ,
#    485 ,
#    1140 ,
#    1191 ,
#    408 ,
#    15 ,
#    354 ,
#    -13 ,
#    423 ]

    #run 920332
    vec_magset = [0,
      -15,
      520,
        0,
      485,
     1140,
     1191,
      408,
       15,
      354,
      -13,
      423,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
454.70001,
    -54.5,
-125.6999,
51.200000,
411.79998,
456.79998,
-0.200000,
1161.5999,
315.29998,
   1543.5,
        0,
        0]

    
    
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


    magset["QFQ1"] = vec_magset[19]
    magset["BFV1"] = vec_magset[20]
    magset["BFH1"] = vec_magset[21]
    magset["BFV2"] = vec_magset[22]
    magset["QFQ2"] = vec_magset[23]
    magset["QFQ3"] = vec_magset[24]
    magset["BFH2"] = vec_magset[25]
    magset["BFVD1"] = vec_magset[26]
    magset["QFQ4"] = vec_magset[27]
    magset["BFVD2"] = vec_magset[28]



    magset["BAD1"] = 0.1 #dummy values for a fake magnet
    magset["BAF1"] = 0.1
    magset["BAD2"] = 0.1
    magset["BAF2"] = 0.1
    magset["BAD3"] = 0.1
    magset["BAF3"] = 0.1
    magset["BAD4"] = 0.1
    magset["BAF4"] = 0.1
    magset["BAD5"] = 0.1
    magset["BAF5"] = 0.1
    magset["BAD6"] = 0.1
    magset["BAF6"] = 0.1
    magset["BAD7"] = 0.1
    magset["BAF7"] = 0.1
    magset["BAD8"] = 0.1
    magset["BAF8"] = 0.1
    magset["BAD9"] = 0.1
    magset["BAF9"] = 0.1
    magset["BAD10"] = 0.1
    magset["BAF10"] = 0.1
    magset["BAD11"] = 0.1
    magset["BAF11"] = 0.1
    magset["BAD12"] = 0.1
    magset["BAF12"] = 0.1
    magset["BAD13"] = 0.1
    magset["BAF13"] = 0.1
    magset["BAD14"] = 0.1
    magset["BAF14"] = 0.1

#    magset["QFQ1"] = 0.1 
#    magset["BFV1"] = 0.1 
#    magset["BFH1"] = 0.1 
#    magset["BFV2"] = 0.1 
#    magset["QFQ2"] = 0.1 
#    magset["QFQ3"] = 0.1 
#    magset["BFH2"] = 0.1 
#    magset["BFVD1"] = 0.1 
#    magset["QFQ4"] = 0.1 
#    magset["BFVD2"] = 0.1 
    kvals = {}

    for magnet in magset:
        if(magnet.startswith("BAF")): #the bending magnets are always set to the same current so can just dfirectly set these here
            #(they need to be set the same each time otherwise the beam will miss the final focus section)
            kvals[magnet] = [-1.61562067, -0.3617/(0.001*spreadsheet.line.at[magnet, 'polelength'])] #B field then K1
            continue
        if(magnet.startswith("BAD")): #the bending angle is slightly different for the focus and defocus, try setting them the the correct B field for their angle
            kvals[magnet] = [-1.528319775, 0.3617/(0.001*spreadsheet.line.at[magnet, 'polelength'])] #B field then K1
            continue
        mag_df = magnet_response[magnet_response['element'] == magnet]
        kvals[magnet] = np.interp(magset[magnet], mag_df['current'], mag_df['kval'])
        zero_field = np.interp(0, mag_df['current'], mag_df['kval'])
        if magnet[0] == 'B': #bending magnets
            kvals[magnet] = -(kvals[magnet]-zero_field) * (proton_momentum/0.2998)  / (0.001*spreadsheet.line.at[magnet, 'polelength'])
        else:
            print(f"{magnet} with field {(kvals[magnet]-zero_field)}")
            kvals[magnet] = (kvals[magnet]-zero_field) / (0.001*spreadsheet.line.at[magnet, 'polelength'])
        if abs(kvals[magnet]) < 1e-3: #if the strength is zero bdsim will treat it as a drift so force it to be non-zero, if its too small the integrator will fall over however
            kvals[magnet] = 1e-3

    if(use_previous_best_fit):
        #overwrite the strengths we fit (allows for partial fitting of the beamline but still including the results in the whole beamline)
        #use a previous fit as the parameters
        if(len(sys.argv) > 1):
            file = ROOT.TFile(sys.argv[1], "READ")
        else:
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




    print(kvals)

    printer = BeamlinePrinter(surv, kvals)

    #exit(0)
    if(use_previous_best_fit):
        printer.print("optimised.gmad")
    else:
        printer.print("unoptimised.gmad")
   

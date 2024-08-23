import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import re
import copy as cp
matplotlib.use('qtagg')

def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

class ssem:
    point1 = []
    point2 = []
    def __init__(self, points):
        if(isinstance(points[0], pd.DataFrame)):
            self.point1 = np.array([points[0]['x2022'].iloc[0], points[0]['y2022'].iloc[0], points[0]['h2022'].iloc[0]])
            self.point2 = np.array([points[1]['x2022'].iloc[0], points[1]['y2022'].iloc[0], points[1]['h2022'].iloc[0]])
        elif(isinstance(points[0], list)):
            self.point1 = cp.deepcopy(np.array(points[0]))
            self.point2 = cp.deepcopy(np.array(points[1]))
        elif(isinstance(points[0], np.ndarray)):
            self.point1 = cp.deepcopy(points[0])
            self.point2 = cp.deepcopy(points[1])
        elif(isinstance(points, ssem)):
            self.point1 = cp.deepcopy(points.point1)
            self.point2 = cp.deepcopy(points.point2)
    def __sub__(self, other):
        tmpp1 = self.point1 - other.point1
        tmpp2 = self.point2 - other.point2
        return ssem([tmpp1, tmpp2])
    def __repr__(self):
        return f"class=ssem, point1='{self.point1}', point2='{self.point2}'"
    def get_center(self):
        center = self.point1
        center[3] = center[3]-225.5
        return center
    def project_along_beamline(self, vec):
        return self.point1.dot(vec) + ssem1s
    def rotatexy(self, angle):
        costh = np.cos(angle)
        sinth = np.sin(angle)
        R = np.array([[costh, sinth, 0],
                    [-sinth, costh, 0],
                    [0, 0, 1]])
        self.point1 = R.dot(self.point1)
        self.point2 = R.dot(self.point2)




def get_ssem(surv, ssemid):
    name1 = 'SSEM'+str(ssemid)+'1'
    name2 = 'SSEM'+str(ssemid)+'2'
    return ssem([surv[surv['name'] == name1], surv[surv['name'] == name2]])

def curvilinear_coords(s):
    start_of_bpd1 = 17400
    intersection_point = start_of_bpd1+3600.5853
    end_of_bpd2 = start_of_bpd1 + 3000 + 1198 + 3000

    if(s<start_of_bpd1):
        x=s
        return x, 0
    elif(s<end_of_bpd2):
        return 0, 0
    else:
        s -= end_of_bpd2
        x = 0.9977756*s + 7193.162 + start_of_bpd1
        y = 0.0666625*s + 240.1859
        return x, y

def coords_after_bpd(mu):
    start_of_bpd1 = 17400
    intersection_point = start_of_bpd1+3600.5853/2.
    P = np.array([intersection_point, 0])
    v = np.array([0.9977756, 0.0666625])
    pos = P + mu*v
    return pos[0], pos[1]


def read_excel(filename):
    surv = pd.ExcelFile(filename)

    surv = surv.parse('基準点成果表 (側壁fit) ', skiprows=3)#.dropna(axis=0, how='all').dropna(axis=1, how='all')
    surv.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters
    surv.columns = ['index', 'ID', 'name', 'x2022', 'y2022', 'h2022', 'x2017', 'y2017', 'h2017', 'x2014', 'y2014', 'h2014', 'dx', 'dy', 'dh', 'comment']
    surv = surv.drop(['x2014', 'y2014', 'h2014'], axis=1)
    surv = surv.drop(['x2017', 'y2017', 'h2017'], axis=1)
    surv = surv.drop(['dx', 'dy', 'dh'], axis=1)
    return surv

def parse_ssem(surv):
    #for now only keep SSEM data
    surv = surv[surv.name.str.contains('SSEM*', regex=True)]
    offset = [[0, 225.5, 0.],[330, 225.5, 0.]] #offset in x,y,s beamline coords

    #SSEM1/2 are God
    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    #coords relative to ssem1[0]
    surv.loc[:,'x2022'] = surv['x2022']-ssem1[0]['x2022'].iloc[0]
    surv.loc[:,'y2022'] = surv['y2022']-ssem1[0]['y2022'].iloc[0]
    surv.loc[:,'h2022'] = surv['h2022']-ssem1[0]['h2022'].iloc[0]

    #estimate initial beamline direction taking ssem2[0]-ssem1[0]
    ssems = [get_ssem(surv, id) for id in range(1,10)]
    #NOTE SSEMs now 0 indexed!!!
    

    delta_ssem10 = ssems[1] - ssems[0]
    print(delta_ssem10.point1)
    beamline_dir = np.array(delta_ssem10.point1)
    print("length along ssem 1->2 ", np.linalg.norm(beamline_dir))
    beamline_dir /= np.linalg.norm(beamline_dir)
    #rotate so that this direction is along the x axis makes calculating curvilinear position easier
    costh = beamline_dir[0]/np.linalg.norm(np.array([beamline_dir[0], beamline_dir[1]])) #normalise to the xy plane
    sinth = beamline_dir[1]/np.linalg.norm(np.array([beamline_dir[0], beamline_dir[1]]))

    angle = -np.arccos(costh)
    [ss.rotatexy(angle) for ss in ssems]


    beamline_dir = np.array([1., 0, beamline_dir[2]])
    beamline_dir /= np.linalg.norm(beamline_dir)
   
    print(beamline_dir)
    print(ssems)
    for i in range(0,9):
        print("SSEM", i+1, "Distance from ssem1 ",ssems[i].project_along_beamline(beamline_dir))
    
    
    return ssems
   # matching_indices = np.argwhere(bi.values == '電流[A]')
   # row = matching_indices[0][0]
   # col = matching_indices[0][1]

   # response = bi.iloc[row:,col:].dropna(axis=0, how='all')
   # response.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters

   # response.columns = response.iloc[0]
   # response = response.drop(response.index[0])
   # response = response[response['[A]']!=0]
   # response['Bmag'] = (response['Bx[Gauss]']**2 + response['By[Gauss]']**2 + response['Bz[Gauss]']**2)**0.5
   # responsemirror = pd.DataFrame({'[A]': -response['[A]'], 'Bmag':  -response['Bmag'], 'Bx[Gauss]': -response['Bx[Gauss]'], 'By[Gauss]': -response['By[Gauss]'], 'Bz[Gauss]': -response['Bz[Gauss]']})
   # response = response.merge(responsemirror, how='outer')
   # return response



pd.set_option('display.max_rows', 500000)

ssem1s = 1769.-85.
survey = read_excel("03_2022_Neutrino.xlsx")
ssems = parse_ssem(survey)

print(ssems[0])

s = np.linspace(0, 50000, 10000)
vec_curvilinear_coords = np.vectorize(curvilinear_coords)
x, y = vec_curvilinear_coords(s)

mu = np.linspace(0, 50000, 1000)
vec_coords_after_bpd = np.vectorize(coords_after_bpd)

x1, y1 = vec_coords_after_bpd(mu)


plt.scatter([s.point1[0] for s in ssems], [s.point1[1] for s in ssems])
plt.scatter(x, y, s=0.2)
plt.scatter(x1, y1, s=0.2)

plt.show()



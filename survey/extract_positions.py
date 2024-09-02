import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import re
import copy as cp
matplotlib.use('qtagg')


bpd1_angle = 0.033547
bpd2_angle = 0.033165

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

def normvec(a):
    retval = np.array(a)
    retval /= np.linalg.norm(retval)
    return retval

def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

class holder:
    point = []
    center = []
    offset = []
    set_center = False
    s = []
    longitudinal = False
    def __init__(self, points, longitudinal=False, offsets=None, center=None, s=np.array([1., 0., 0.])):
        self.s = s.T/np.linalg.norm(s)
        self.longitudinal = longitudinal
        if(isinstance(points[0], pd.DataFrame)):
            self.point = [np.array([pt['x2022'].iloc[0], pt['y2022'].iloc[0], pt['h2022'].iloc[0]]) for pt in points]
        elif(isinstance(points[0], list)):
            self.point = [cp.deepcopy(np.array(pt)) for pt in points]
        elif(isinstance(points[0], np.ndarray)):
            self.point = [cp.deepcopy(pt) for pt in points]
        elif(isinstance(points, holder)):
            self.point = cp.deepcopy(points.point)
        if center is not None:
            self.center = center
            self.set_center = True
        if offsets is not None:
            self.offset = cp.deepcopy(offsets)
            self.calculate_center()
        else:
            self.offset = np.nan

    def __sub__(self, other):
        return holder([self.point[i] - other.point[i] for i in range(len(self.point))], center = self.get_center() - other.get_center())
    def __repr__(self):
        return f"class=holder, point='{self.point}', center='{self.center}'"

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

    def calculate_center(self):  ##TODO make this work demanding that the point[0]->center vector has no component along the beamline...
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
#        phi = np.linspace(0, 2.*np.pi, 100)
#        circ = np.array([locus(invtanopen, R2, conelen, ang) for ang in phi])


#        #now rotate that vector about vec until the resulting vector has the correct component along the beamline
#        costh = p1c.dot(-vec)
#        sinth = (1-costh**2)**0.5
#        invtanth = costh/sinth
#        print(invtanth)
#        #now get rotation matrix for a 001 vector to direction -vec
#        openangle = invtanth * np.ones(100)
#        cone = np.array([np.cos(phi), np.sin(phi), openangle]).T
#        cone = np.array([c/np.linalg.norm(c) for c in cone])
#        rotcone = R2.dot(cone.T).T
#        centerline_to_cone = self.point[0] - (self.point[1] + np.linalg.norm(self.point[1] - self.point[0]) * rotcone)
#
#
#        center_points = self.point[0] + np.linalg.norm(self.offset[0]) * rotcone
#

        vertical = np.array([0., 0., 1.])
        s = self.s
        if(self.longitudinal):
            s = np.cross(vertical, self.s)


#        proj = lambda invtan, R, s, theta: self.point[0] - (self.point[1] + np.linalg.norm(self.point[1] - self.point[0]) * R.dot(normvec([np.cos(theta), np.sin(theta), invtan]).T))
#
#
#        projection = centerline_to_cone.dot(self.s)
#        newproj = [proj(invtanth, R2, self.s, ang) for ang in phi]
#        print(np.linalg.norm(np.array(centerline_to_cone[0])))
#        print(np.linalg.norm(newproj[0]))
#        exit(1)
#
#        plt.plot(phi, projection)
#        plt.plot(phi, newproj)
#        plt.show()
        

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

def get_holder(surv, basename, points, longitudinal, offsets, s=None):
    names = [basename+str(pt) for pt in points]
    if(s is None):
        return holder([surv[surv['name'] == name] for name in names],  longitudinal, offsets)
    else:
        return holder([surv[surv['name'] == name] for name in names],  longitudinal, offsets, s)

def curvilinear_coords(s):
    start_of_bpd1 = 17400.
    intersection_point = start_of_bpd1+3600.5853
    end_of_bpd2 = start_of_bpd1 + 3000. + 1198. + 3000.

    if(s<start_of_bpd1):
        x=s
        return x, 0.
    elif(s<end_of_bpd2):
        return 0., 0.
    else:
        s -= end_of_bpd2
        x = np.cos(bpd1_angle+bpd2_angle)*s + 7193.162 + start_of_bpd1
        y = np.sin(bpd1_angle+bpd2_angle)*s + 240.1859
        return x, y

def coords_after_bpd(mu):
    start_of_bpd1 = 17400
    intersection_point = start_of_bpd1+3600.5853
    P = np.array([intersection_point, 0])
    v = np.array([0.9977756, 0.0666625])
    pos = P + mu*v
    return pos[0], pos[1]


def read_excel(filename):
    surv = pd.ExcelFile(filename)
    surv = surv.parse('基準点成果表 (側壁fit) ', skiprows=3)
    surv.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters
    surv.columns = ['index', 'ID', 'name', 'x2022', 'y2022', 'h2022', 'x2017', 'y2017', 'h2017', 'x2014', 'y2014', 'h2014', 'dx', 'dy', 'dh', 'comment']
    surv = surv.drop(['x2014', 'y2014', 'h2014'], axis=1)
    surv = surv.drop(['x2017', 'y2017', 'h2017'], axis=1)
    surv = surv.drop(['dx', 'dy', 'dh'], axis=1)
    surv.loc[surv['name'] == 'ssem12', 'h2022'] += 50. #CERN TODO REMOVE!!!
    return surv

def rotate_surveyxy(surv, R):
    tmp = surv
    surv['x2022'] = R[0,0] * tmp['x2022'] + R[0,1] * tmp['y2022']
    surv['y2022'] = R[1,0] * tmp['x2022'] + R[1,1] * tmp['y2022']
    return surv

def parse_holder(surv, ssem1s):
    #for now only keep SSEM data
#    surv = surv[surv.name.str.contains('SSEM*', regex=True)]
    ssem_offset = np.array([[0., 0., 225.5],[0., 330., 225.5]]) #offset in x,y,h where x is along the beamline,y is horisontal and h is vertical
    bpd_offset =  np.array([[0., 0., 500.],[3000., 0., 500.]]) #offset in x,y,h where x is along the beamline,y is horisontal and h is vertical
    qpq_offset =  np.array([[0., 0., 500.],[3000., 0., 500.]]) #offset in x,y,h where x is along the beamline,y is horisontal and h is vertical

    #due to the issue of multiple solutions, described in the rotation matrix we must make an assumption:
    #components are placed such that the vertical axis of their vertical face is perpendicular to the beamline 
    #essentially the points are on top of the components, so allow them to hang from these points and be pulled vertical by gravity


    #SSEM1/2 are God
    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    ssem2 = [surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']]
   
    #estimate initial beamline direction taking ssem2[0]-ssem1[0]
    s = np.array([ssem2[0][i].iloc[0] - ssem1[0][i].iloc[0] for i in ['x2022', 'y2022', 'h2022']])

    ssem1_center = holder(ssem1, False, ssem_offset, s=s).get_center()

    #coords relative to ssem1[0]
    surv.loc[:,'x2022'] = surv['x2022']-ssem1_center[0]
    surv.loc[:,'y2022'] = surv['y2022']-ssem1_center[1]
    surv.loc[:,'h2022'] = surv['h2022']-ssem1_center[2]


    theta = np.atan(s[1]/s[0])
    theta += np.pi #we know that the survey is in the bottom-right quadrant
    R = np.array([[np.cos(theta), -np.sin(theta)],
                [np.sin(theta), np.cos(theta)]])
    surv = rotate_surveyxy(surv, R)

    ssem1 = [surv[surv['name'] == 'SSEM11'], surv[surv['name'] == 'SSEM12']]
    ssem2 = [surv[surv['name'] == 'SSEM21'], surv[surv['name'] == 'SSEM22']]

    s = np.array([ssem2[0][i].iloc[0] - ssem1[0][i].iloc[0] for i in ['x2022', 'y2022', 'h2022']])

    s = np.array([0.8, 0.2, 0.])

    ssems = [get_holder(surv, 'SSEM'+str(id), [1, 2], False, ssem_offset, s=s) for id in range(1,10)]
    qpqs = [get_holder(surv, 'PQ'+str(id), [1, 2], True, qpq_offset, s=s) for id in range(1, 6)]
    bpds = [get_holder(surv, 'PD'+str(id), [1, 2], True, bpd_offset, s=s) for id in range(1, 3)]

    #NOTE SSEMs now 0 indexed!!!

    print("ssem1 position ", ssems[0].get_center())
    print("ssem1 points ", ssems[0])
    beamline_dir = ssems[1].get_center() - ssems[0].get_center()
    print("length along ssem 1->2 ", np.linalg.norm(beamline_dir))
    beamline_dir /= np.linalg.norm(beamline_dir)
    R = get_rotation_matrix(beamline_dir, np.array([1., 0., 0.]))
    #rotate so that this direction is along the x axis makes calculating curvilinear position easier


    print("qpq1 before rotation ", qpqs[0])
#    [ss.rotate(R) for ss in ssems] TODO uncomment
#    [bp.rotate(R) for bp in bpds]
#    [qp.rotate(R) for qp in qpqs]
#    print("qpq1 after rotation ", qpqs[0])
#    offset = np.array([ssem1s, 0., 0.])
#    [ss.set_zero(offset) for ss in ssems]
#    [bp.set_zero(offset) for bp in bpds]
#    [qp.set_zero(offset) for qp in qpqs]
 
    print("qpq1 after set zero ", qpqs[0])

    print("ssem1 position after rotation ", ssems[0].get_center())
    print("ssem1 point1 after rotation ", ssems[0].point[0])

    #beamline_dir = np.array([1., 0, beamline_dir[2]])
    #beamline_dir /= np.linalg.norm(beamline_dir) 

    
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
        ax.scatter(x, y, z, label=lab)
    else:
        plt.scatter(x, y, label=lab)
      

pd.set_option('display.max_rows', 500000)

ssem1s = 1769.-85.
survey = read_excel("03_2022_Neutrino.xlsx")
ssems, bpds, qpqs = parse_holder(survey, ssem1s)

s = np.linspace(0, 50000, 10000)
vec_curvilinear_coords = np.vectorize(curvilinear_coords)
x, y = vec_curvilinear_coords(s)

mu = np.linspace(0, 50000, 1000)
vec_coords_after_bpd = np.vectorize(coords_after_bpd)

x1, y1 = vec_coords_after_bpd(mu)

print("qpq1 before plot1 ", qpqs[0])

plot_points(ssems, 'SSEM', 0, 2)
plot_points(bpds, 'BPD', 0, 2)
plot_points(qpqs, 'QPQ', 0, 2)
plt.xlabel('x(mm)')
plt.ylabel('z(mm)')
plt.show()


plot_centers(ssems, 'SSEM', 0, 2)
plot_centers(bpds, 'BPD', 0, 2)
plot_centers(qpqs, 'QPQ', 0, 2)
plt.xlabel('x(mm)')
plt.ylabel('z(mm)')
plt.show()




plot_points(ssems, 'SSEM', 0, 1)
plot_points(bpds, 'BPD', 0, 1)
plot_points(qpqs, 'QPQ', 0, 1)

plot_centers(ssems, 'SSEM center', 0, 1)
plot_centers(bpds, 'BPD center', 0, 1)

print("qpq1 before plot2 ", qpqs[0])
plot_centers(qpqs, 'QPQ center', 0, 1)

#plt.scatter([s.point[0][0] for s in ssems], [s.point[0][1] for s in ssems])
#plt.scatter([b.point[0][0] for b in bpds], [b.point[0][1] for b in bpds])
#plt.scatter([b.point[1][0] for b in bpds], [b.point[1][1] for b in bpds])


plt.scatter(x, y, s=0.2)
plt.scatter(x1, y1, s=0.2)
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




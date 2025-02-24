import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import re
import copy as cp
import sys
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
misalignments=False
bias_physics=False
print_vacuum=False


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


def string_to_list(s):
    s = s.strip('[]').strip()
    elements = s.split()
    return elements


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
        if(isinstance(points[0], pd.DataFrame)):
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
            self.offset = cp.deepcopy(offsets)
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
        
        self.ssems, self.bpds, self.qpqs = parse_holder(self.survey, ssem1s)
        self.line['s_start'] = self.line['length'].shift().cumsum()

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
            if not np.isnan(row['s_survey']) and row['type'] == 'ssem':
                x.append(row.s_start)
                y.append(row.s_start-row.s_survey)
                print(row['element'],',\t\t\t',"{:.3f}".format(row['s_start']+85.), ',\t',"{:.3f}".format(row['s_survey']), ',\t',"{:.3f}".format(row['s_start']-row['s_survey']+85.0))
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
        if(not geometry):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.print_xsec_bias('vacuum')
        self.endl()

    def print_quad_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, tilt='+str(row.tilt)+', k1='+str(self.kvals[row.element]))
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
            
            if(row.type != 'drift' and previously_drift and not geometry):
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
                if(geometry):
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
        return holder([surv[surv['name'] == name] for name in names], longitudinal, offsets)
    else:
        return holder([surv[surv['name'] == name] for name in names],  longitudinal, offsets, s)

def curvilinear_coords(s):
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
    qpq_svec = [s if i<=2 else s_after_bpd for i in range(1,6)]

    ssems = [get_holder(surv, 'SSEM'+str(id), [1, 2], False, ssem_offset, s=ssem_svec[id-1]) for id in range(1,10)]
    qpqs = [get_holder(surv, 'PQ'+str(id), [1, 2], True, qpq_offset, s=qpq_svec[id-1]) for id in range(1, 6)]
    bpds = [get_holder(surv, 'PD'+str(id), [1, 2], True, bpd_offset, s=s) for id in range(1, 3)]

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


    exit(1)


    survey = read_excel("03_2022_Neutrino.xlsx")
    ssems, bpds, qpqs = parse_holder(survey, ssem1s)


    s = np.linspace(0, 50000, 10000)
    vec_curvilinear_coords = np.vectorize(curvilinear_coords)
    x, y = vec_curvilinear_coords(s)

    mu = np.linspace(0, 50000, 1000)


    plot_points(ssems, 'SSEM Survey', 0, 2)
    plot_points(bpds, 'BPD Survey', 0, 2)
    plot_points(qpqs, 'QPQ Survey', 0, 2)
    plot_centers(ssems, 'SSEM center', 0, 2)
#    plot_centers(bpds, 'BPD center', 0, 2)
#    plot_centers(qpqs, 'QPQ center', 0, 2)
    plt.xlabel('x(mm)')
    plt.ylabel('z(mm)')
    plt.legend()
    plt.show()


    plot_points(ssems, 'SSEM', 0, 1)
    plot_points(bpds, 'BPD', 0, 1)
    plot_points(qpqs, 'QPQ', 0, 1)
    plot_centers(ssems, 'SSEM center', 0, 1)
#    plot_centers(bpds, 'BPD center', 0, 1)
#    plot_centers(qpqs, 'QPQ center', 0, 1)
    plt.scatter(x, y, s=0.2)
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
#    plot_centers(bpds, 'BPD center', 0, 1, 2, ax)
    plot_points(qpqs, 'QPQ', 0, 1, 2, ax)
#    plot_centers(qpqs, 'QPQ center', 0, 1, 2, ax)
    ax.set_xlabel('x(mm)')
    ax.set_ylabel('y(mm)')
    ax.set_zlabel('z(mm)')

    plt.legend()
    #ax.set_aspect('equal', 'box')
    plt.show()




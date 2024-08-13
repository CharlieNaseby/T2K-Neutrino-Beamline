import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

def string_to_list(s):
    s = s.strip('[]').strip()
    elements = s.split()
    return elements

def extract_number(string):
    # Define the regex pattern to match a string with a number in the middle
    pattern = r'\D*(\d+)\D*'
    # Search for the pattern in the string
    match = re.search(pattern, string)
    # If a match is found, return the number, otherwise return None
    if match:
        return int(match.group(1))
    else:
        return None



class BeamlinePrinter:
    def __init__(self, line, filename, primaries_only=False):
        self.beamline = line
        self.file = open(filename, "w")
        self.s = 0
        self.blmID = 1
        self.primaries_only=primaries_only
        self.line = []

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
        self.file.write(';\n')

    def print_bend_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, angle='+str(row.angle)+', tilt='+str(row.tilt)+', B='+str(kvals[row.element])+'*T')
        if(no_geom):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.file.write(';\n')

    def print_quad_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, tilt='+str(row.tilt)+', k1='+str(kvals[row.element]))
        if(no_geom):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.file.write(';\n')

    def print_field(self, row, ndim):
         self.file.write(row.element+'field: field, type="bmap'+str(ndim)+'d", bScaling=1.0, magneticFile="bdsim'+str(ndim)+'d:../magnet_responses/'+row.element+'.dat";\n')
       
    def print_fieldmapgeom(self, row, ndim):
        self.line.append(row.element)
        self.print_field(row, ndim)
        self.file.write(row.element+': element, geometryFile="gdml:../'+row.element+'.gdml", fieldAll="'+row.element+'field", l='+str(row.length)+'*mm;\n')

    def print_fieldmap(self, row, magtype, ndim):
        self.line.append(row.element)
        self.print_field(row, ndim)
        self.file.write(row.element+': '+magtype+', fieldVacuum="'+row.element+'field", l='+str(row.length)+'*mm, angle='+str(row.angle))
        if(no_geom):
            self.file.write(', magnetGeometryType="none"')
        self.print_aperture(row)
        self.file.write(';\n')

    def print_blm(self, reference, dx, dy, ds, orientation):
        self.file.write('blm_'+reference+'_'+str(self.blmID)+': blm, scoreQuantity="chrg eDep", geometryType="cylinder", blm1=100*mm, blm2=30*mm, blmMaterial="Al",')
        self.file.write('referenceElement="'+reference+'", x='+str(dx)+'*mm, y='+str(dy)+'*mm, s='+str(ds)+'*mm')
        if(orientation=='perp'):
            self.file.write(', theta=1.570796, psi=1.570796;\n')
        else:
            self.file.write(';\n')
        self.blmID+=1

    def print_blms(self, row):
        dx = string_to_list(row.blm_offset_x)
        dy = string_to_list(row.blm_offset_y)
        ds = [str(float(entry)+row.polelength/2.0) for entry in string_to_list(row.blm_offset_s)]  #TODO THIS ONLY WORKS FOR PURE QUAD/DIPOLE FIELDMAP WILL BREAK THIS
        orientation = string_to_list(row.blm_orientation)
        for i in range(len(dx)):
            self.print_blm(row.element, dx[i], dy[i], ds[i], orientation[i])

    def print_ssem(self, row, thickness):
        first_driftlen = row.mark
        name = row.element + "_udrift"
        self.print_drift(row, name, row.mark)
        self.print_collimator(row.element, str(thickness)+"*mm", "G4_Ti", row.aperture_x)
        name = row.element + "_ddrift"
        self.print_drift(row, name, row.length - float(row.mark) - thickness) 

    def print_collimator(self, name, length, material, hWidth):
        self.line.append(name)
        self.file.write(name+": target, l="+str(length)+", material=\""+material+"\", horizontalWidth="+str(hWidth)+"*mm;\n")

    def print_dump(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': dump, horizontalWidth='+str(row.aperture_x)+'*mm, l='+str(row.length)+'*mm')
        self.file.write(';\n')

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

    def print_physics(self, physics_list):
        self.file.write('option, physicsList =' + '"' + physics_list + '";\n')
    
    def sample_ssems(self):
        for element in self.line:
            if(re.match('SSEM[0-9]$', element)):
                self.file.write('sample, range='+element+';\n')

    def print(self):
        self.file.write('chrg: scorer, type="cellcharge";\n')
        self.file.write('eDep: scorer, type="depositeddose";\n')
        for row in self.beamline.itertuples():

            if(row.type in ['rbend', 'sbend', 'quadrupole']):
                self.file.write('\n') #give some breathing room
                self.print_drift(row, row.element+'_driftu', float(row.mark)-row.polelength/2.)
                if(row.type == 'quadrupole'):
                    self.print_quad_magnet(row)
                else:
                    self.print_bend_magnet(row)
                self.print_drift(row, row.element+'_driftd', row.length - (float(row.mark)+row.polelength/2.))
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
            else:
                self.print_drift(row, row.element, row.length)
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
            self.print_beam_from_file("../run_0910580_100k.root")
        elif(beam_halo):
            self.print_halo()
        else:
            self.print_beam_0910580()

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


proton_momentum = 30.924 # momentum for a 30GeV KE proton 

#run 910580
vec_magset = [0 ,
-15 ,
520 ,
0 ,
485 ,
1139 ,
1190 ,
408 ,
15 ,
354 ,
-13 ,
423 ]
  

beamline = strip_whitespace(pd.read_csv("fujii-san.csv", header=0, skipinitialspace=True))
magnet_response = strip_whitespace(pd.read_csv("kicurve.csv", header=0, skipinitialspace=True))


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
        kvals[magnet] = -(kvals[magnet]-zero_field) * (proton_momentum/0.2998)  / (0.001*beamline.loc[beamline['element'] == magnet].iloc[0]['polelength'])
    else:
        kvals[magnet] = (kvals[magnet]-zero_field) / (0.001*beamline.loc[beamline['element'] == magnet].iloc[0]['polelength'])


print(kvals)
#kvals['BPD1'] = -1.15329
#kvals['BPD2'] = -1.14018
#kvals['QPQ4'] = -0.0518735 #set the QPQ4 val to that in the fake fieldmap

print(kvals)
mag_df = magnet_response[magnet_response['element'] == 'BPH3']
#plt.scatter(mag_df['current'], mag_df['kval'])
#plt.show()



print_tunnel=False
print_physics=False
sample_all=False
sample_ssem=True
sample_entry=False
beam_from_file = True
beam_halo = False
enable_blms = False
no_geom = True

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



#df_tmp = magnet_response[magnet_response['element'] == row.element]


prnt = BeamlinePrinter(beamline, "gmad/test.gmad", primaries_only=generate_primaries)
prnt.print()
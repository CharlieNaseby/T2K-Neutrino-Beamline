import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df


class BeamlinePrinter:
    def __init__(self, line, filename):
        self.beamline = line
        self.file = open(filename, "w")
        self.s = 0
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
        self.print_aperture(row)
        self.file.write(';\n')

    def print_quad_magnet(self, row):
        self.line.append(row.element)
        self.file.write(row.element+': '+row.type+', l='+str(row.polelength)+'*mm, tilt='+str(row.tilt)+', k1='+str(kvals[row.element]))
        self.print_aperture(row)
        self.file.write(';\n')
    
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
        self.file.write(row.element+': dump, horizontalWidth='+str(row.aperture_x)+'*mm')
        self.file.write(';\n')

    def print_tunnel(self):
        self.file.write('''
option, buildTunnel = 1,
tunnelType="rectangular",
tunnelOffsetX = 100*cm,
tunnelOffsetY = 100*cm,
tunnelAper1 = 150*cm,
tunnelAper2 = 150*cm,
tunnelThickness = 30*cm,
buildTunnelFloor = 0,
tunnelSoilThickness = 2*m;\n\n''')

    def print_physics(self, physics_list):
        self.file.write('option, physicsList =' + '"' + physics_list + '";\n')


    def print(self):
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
            elif(row.type == 'ssem'):
                self.print_ssem(row, 15e-3)
            elif(row.type == 'wsem'):
                self.print_ssem(row, 1e-6)
            elif(row.type == 'dump'):
                self.print_dump(row)
            else:
                self.print_drift(row, row.element, row.length)
            self.s += row.length
            self.file.write("! s=" + str(self.s) + "\n")

        #print the beamline elements in a line
        self.file.write('\nl0: line = (')
        for idx, element in enumerate(self.line):
            if idx == len(self.line)-1:
                self.file.write(element +');\n')
            else:
                self.file.write(element +',\n')

        self.file.write('\nuse, period=l0;\n')
        self.print_beam()
        self.print_tunnel()
        self.print_physics('g4FTFP_BERT')
        self.file.write('option, nturns=1;')
        self.file.write('sample, all;\n')


#bfield = {'BPV1': 0.0001,
#          'BPH2': 0.0001,
#          'BPD1': 1.15329,
#          'BPD2': 1.14018,
#          'BPV2': 0.0001,
#          'BPH3': 0.0001,}
#
#k1 = {'QPQ1': 0.0511454667,
#      'QPQ2': -0.060356033,
#      'QPQ3': 0.07795667,
#      'QPQ4': -0.0132418,
#      'QPQ5': 0.09191465}

proton_momentum = 30.924 # momentum for a 30GeV KE proton 

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
8 ,
423]


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


kvals['BPD1'] = -1.15329
kvals['BPD2'] = -1.14018

print(kvals)
mag_df = magnet_response[magnet_response['element'] == 'BPH3']
#plt.scatter(mag_df['current'], mag_df['kval'])
#plt.show()


#df_tmp = magnet_response[magnet_response['element'] == row.element]


prnt = BeamlinePrinter(beamline, "gmad/test.gmad")
prnt.print()
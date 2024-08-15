import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

def include_sign(df, x, b):
    zeroval = df[df[x] == 0]
    df[b] = df[b]*np.sign(df[x])
    dftmp = df
    dftmp['mag'] = abs(dftmp[x])
    closest = dftmp.sort_values(by='mag').head(3)
    neg = closest[closest[x] < 0.]
    pos = closest[closest[x] > 0.]
    #average these points to get exp value at 0
    zero_sign = np.sign(neg[b].iloc[0] + pos[b].iloc[0])
    df.loc[df[x] == 0, b] = zeroval[b].iloc[0]*zero_sign
    return df

def interpolate_point(df, val, incol, outcol):
    points = df.iloc[(df[incol]-val).abs().argsort()[:2]] 
    x=np.array(points.iloc[:][incol])
    if(val > x.max() or val < x.min()):
        return points.iloc[0][outcol]
    else:
        return(points.iloc[0][outcol] + (points.iloc[1][outcol]-points.iloc[0][outcol])/(points.iloc[1][incol] - points.iloc[0][incol]) * (val - points.iloc[0][incol]))
   
def interpolate(df, val, incol, outcol):
    return [interpolate_point(df, v, incol, outcol) for v in val]

def read_bi(filename):
    bi = pd.ExcelFile(filename)

    bi = bi.parse(0).dropna(axis=0, how='all').dropna(axis=1, how='all')

    matching_indices = np.argwhere(bi.values == '電流[A]')
    row = matching_indices[0][0]
    col = matching_indices[0][1]

    response = bi.iloc[row:,col:].dropna(axis=0, how='all')
    response.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters

    response.columns = response.iloc[0]
    response = response.drop(response.index[0])
    response = response[response['[A]']!=0]
    response['Bmag'] = (response['Bx[Gauss]']**2 + response['By[Gauss]']**2 + response['Bz[Gauss]']**2)**0.5
    responsemirror = pd.DataFrame({'[A]': -response['[A]'], 'Bmag':  -response['Bmag'], 'Bx[Gauss]': -response['Bx[Gauss]'], 'By[Gauss]': -response['By[Gauss]'], 'Bz[Gauss]': -response['Bz[Gauss]']})
    response = response.merge(responsemirror, how='outer')
    return response

def read_xyz(infile, colname):
    z = pd.ExcelFile(infile)

    z = z.parse(0)
    z.replace({r'[^\x00-\x7F]+':''}, regex=True, inplace=True) #remove japanese characters
    z.columns = z.iloc[0]
    z = z.drop(z.index[0])
    z1500 = z.iloc[:, 0:2]
    z1300 = z.iloc[:, 2:4]
    z1100 = z.iloc[:, 4:6]

    ratio = z1500[colname] / z1100[colname]
    error = ((0.5/z1500[colname])**2 + (0.5/z1100[colname])**2)**0.5
    z1100 = z1100.dropna()
#    plt.errorbar(z1500['Z[mm]'], ratio, yerr=error, lw=0, elinewidth=2)
#    plt.scatter(z1500['Z[mm]'], ratio)
#    plt.show()
    return z1100

def create_fieldmap3d(magnet, magset, magset_ref, sign):
    print(magnet)
    x = read_xyz("magnet_responses/"+magnet+"/"+magnet+"_X.xls", 'By[Gauss]')
    y = read_xyz("magnet_responses/"+magnet+"/"+magnet+"_Y.xls", 'Bx[Gauss]')
    z = read_xyz("magnet_responses/"+magnet+"/"+magnet+"_Z.xls", 'By[Gauss]')
    bi = read_bi("magnet_responses/"+magnet+"/"+magnet+"_BI.xls")
    
    nzpoints = 20 #how many z positions to include in the half magnet
    length = float(beamline.loc[beamline['element'] == magnet, 'length'].item())
    polelen = float(beamline.loc[beamline['element'] == magnet, 'polelength'].item()) #the length of the pole of the magnet in mm
    mark = float(beamline.loc[beamline['element'] == magnet, 'mark'].item())
    zoffset = mark - 0.5*length
    fringelen = min(length - (mark+polelen/2.), mark-polelen/2.) #length of the fringe field to include outside of the magnet face in mm
    zref = float(beamline.loc[beamline['element'] == magnet, 'zref'].item()) #z position at which the x and y scans were made
    fieldmap_current = magset_ref
    real_current = magset
    
    #get the field strength scaling to apply to the x, y, z fieldmaps
    field_scale = interpolate_point(bi, real_current, '[A]', 'By[Gauss]')/interpolate_point(bi, fieldmap_current, '[A]', 'By[Gauss]')
    x['By[Gauss]'] = field_scale*x['By[Gauss]']
    y['Bx[Gauss]'] = field_scale*y['Bx[Gauss]']
    z['By[Gauss]'] = field_scale*z['By[Gauss]']
    x = include_sign(x, 'X[mm]', 'By[Gauss]')
    y = include_sign(y, 'Y[mm]', 'Bx[Gauss]')
    x['By[Gauss]'] = sign*x['By[Gauss]']
    y['Bx[Gauss]'] = sign*y['Bx[Gauss]']
   
    zfield = pd.DataFrame({'z': np.linspace(-fringelen, polelen/2., nzpoints)})
    zfield['By[Gauss]'] = interpolate(z, zfield['z'], 'Z[mm]', 'By[Gauss]')
    zfieldmirror = pd.DataFrame({'z': polelen-zfield['z'], 'By[Gauss]':  zfield['By[Gauss]']})
    zfield = zfield.merge(zfieldmirror, how='outer')
    
    zscaling = zfield
    zscaling['scale'] = zscaling['By[Gauss]']/interpolate_point(z, zref, 'Z[mm]', 'By[Gauss]')
    zscaling.drop(['By[Gauss]'], axis=1)
    y['Fz'] = 0
    
    #fieldmap is in cm
    x['X[cm]'] = x['X[mm]']/10.
    y['Y[cm]'] = y['Y[mm]']/10.
    zscaling['z'] = (zscaling['z'] - polelen/2. + zoffset)/10.
    
    xy = y.assign(key=1).merge(x.assign(key=1), how='outer', on='key')
    
    outf = open("magnet_responses/"+magnet+".dat", "w")
    
    outf.write('xmin> '+ str(x['X[cm]'].min())+'\n')
    outf.write('xmax> '+ str(x['X[cm]'].max())+'\n')
    outf.write('nx> '+ str(x.shape[0])+'\n')
    
    outf.write('ymin> '+ str(y['Y[cm]'].min())+'\n')
    outf.write('ymax> '+ str(y['Y[cm]'].max())+'\n')
    outf.write('ny> '+ str(y.shape[0])+'\n')
    
    outf.write('zmin> '+ str(zscaling['z'].min())+'\n')
    outf.write('zmax> '+ str(zscaling['z'].max())+'\n')
    outf.write('nz> '+ str(zscaling.shape[0])+'\n')
    
    outf.write('! X \t Y \t Z \t Fx \t Fy \t Fz'+'\n')
    
    for indx, zslice in zscaling.iterrows():
        tmp = xy
        tmp['Fx'] = tmp['Bx[Gauss]']*zslice['scale']/10000. 
        tmp['Fy'] = tmp['By[Gauss]']*zslice['scale']/10000.
        tmp['Fz'] = tmp['Fz']*zslice['scale']/10000.
        tmp['z'] = zslice['z']
        for ind, row in tmp.iterrows():
            outf.write(str(row['X[cm]'])+ " "+ str(row['Y[cm]'])+" "+str(row['z'])+" "+str(row['Fx'])+" "+str(row['Fy'])+" "+str(row['Fz'])+"\n")
    outf.close()


def create_fieldmap1d(magnet, magset, magset_ref, sign, horisontal_bending):
    print(magnet)
    z = read_xyz("magnet_responses/"+magnet+"/"+magnet+"_Z.xls", 'By[Gauss]')
    bi = read_bi("magnet_responses/"+magnet+"/"+magnet+"_BI.xls")

    nzpoints = 20 #how many z positions to include in the half magnet
    length = float(beamline.loc[beamline['element'] == magnet, 'length'].item())
    polelen = float(beamline.loc[beamline['element'] == magnet, 'polelength'].item()) #the length of the pole of the magnet in mm
    mark = float(beamline.loc[beamline['element'] == magnet, 'mark'].item())
    zoffset = mark - 0.5*length
    fringelen = min(length - (mark+polelen/2.), mark-polelen/2.) #length of the fringe field to include outside of the magnet face in mm
    zref = float(beamline.loc[beamline['element'] == magnet, 'zref'].item()) #z position at which the x and y scans were made
    fieldmap_current = magset_ref
    real_current = magset

    #get the field strength scaling to apply to the x, y, z fieldmaps
    print("real field ", interpolate_point(bi, real_current, '[A]', 'By[Gauss]'))
    print("fieldmap field", interpolate_point(bi, fieldmap_current, '[A]', 'By[Gauss]'))
    field_scale = interpolate_point(bi, real_current, '[A]', 'By[Gauss]')/interpolate_point(bi, fieldmap_current, '[A]', 'By[Gauss]')
    z['By[Gauss]'] = field_scale*z['By[Gauss]']


    zfield = pd.DataFrame({'z': np.linspace(-fringelen, polelen/2., nzpoints)})
    zfield['By[Gauss]'] = interpolate(z, zfield['z'], 'Z[mm]', 'By[Gauss]')
    zfieldmirror = pd.DataFrame({'z': polelen-zfield['z'], 'By[Gauss]':  zfield['By[Gauss]']})
    zfield = zfield.merge(zfieldmirror, how='outer')
    
    zscaling = zfield
    if(real_current == 0):
        zscaling['scale'] = 0
    else:
        zscaling['scale'] = zscaling['By[Gauss]']/interpolate_point(z, zref, 'Z[mm]', 'By[Gauss]') #TODO deal with the case that By is zero in both df
    
    #fieldmap is in cm
    zscaling['z'] = (zscaling['z'] - polelen/2. + zoffset)/10.
    
    outf = open("magnet_responses/"+magnet+".dat", "w")
    
    outf.write('zmin> '+ str(zscaling['z'].min())+'\n')
    outf.write('zmax> '+ str(zscaling['z'].max())+'\n')
    outf.write('nz> '+ str(zscaling.shape[0])+'\n')
    
    outf.write('! Z \t Fx \t Fy \t Fz'+'\n')
    
    for indx, zslice in zscaling.iterrows():
        tmp = zslice
        if(horisontal_bending):
            tmp['Fx'] = 0. #tmp['Bx[Gauss]']*zslice['scale']/10000. 
            tmp['Fy'] = sign*tmp['By[Gauss]']*tmp['scale']/10000.
        else:
            tmp['Fy'] = 0. #tmp['Fz']*zslice['scale']/10000.
            tmp['Fx'] = sign*tmp['By[Gauss]']*tmp['scale']/10000.

        tmp['Fz'] = 0. #tmp['Fz']*zslice['scale']/10000.
        tmp['z'] = zslice['z']
        outf.write(str(tmp['z'])+" "+str(tmp['Fx'])+" "+str(tmp['Fy'])+" "+str(tmp['Fz'])+"\n")
    outf.close()


pd.set_option('display.max_rows', 500000)
#vec_magset = [0 ,
#-15 ,
#520 ,
#0 ,
#485 ,
#1139 ,
#1190 ,
#408 ,
#15 ,
#354 ,
#-13 ,
#423 ]

vec_magset = 100*np.ones(12) #set magnet current to 100A

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

beamline = strip_whitespace(pd.read_csv("fujii-san.csv", header=0, skipinitialspace=True))

beamline = beamline[beamline['polelength'] != 0] #get all magnets

quads = ['QPQ1', 'QPQ2', 'QPQ4', 'QPQ5']
dipoles = ['BPD2', 'BPV1', 'BPV2', 'BPH3']

fieldmap_current = {"BPV1": 1800,
                    "BPH2": 1800,
                    "QPQ1": 1012.6,
                    "QPQ2": 1500,
                    "BPD2": 2000,
                    "BPV2": 1500,
                    "BPH3": 1500,
                    "QPQ4": 1100,
                    "QPQ5": 1500}

#are these focusing or defocusing quads
signs = {"QPQ1": -1,
         "QPQ2": 1,
         "QPQ3": -1,
         "QPQ4": 1,
         "QPQ5": -1,
         "BPV1": 1,
         "BPH3": 1,
         "BPD2": -1,
         "BPV2": 1}

horisontal_bending = {"BPV1": False,
                      "BPH2": True,
                      "BPH3": True,
                      "BPD2": True,
                      "BPV2": False}


[create_fieldmap3d(magnet, magset[magnet], fieldmap_current[magnet], signs[magnet]) for magnet in quads] 
[create_fieldmap1d(magnet, magset[magnet], fieldmap_current[magnet], signs[magnet], True) for magnet in dipoles] 


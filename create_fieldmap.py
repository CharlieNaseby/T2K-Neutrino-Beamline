import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def include_sign(df, x, b, sign):
    zeroval = df[df[x] == 0]
    df[b] = df[b]*sign*np.sign(df[x])
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

    response['Bmag'] = (response['Bx[Gauss]']**2 + response['By[Gauss]']**2 + response['Bz[Gauss]']**2)**0.5
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
#    plt.errorbar(z1500['Z[mm]'], ratio, yerr=error, lw=0, elinewidth=2)
#    plt.scatter(z1500['Z[mm]'], ratio)
#    plt.show()
    return z1100


x = read_xyz("magnet_responses/PQ4_X.xls", 'By[Gauss]')
y = read_xyz("magnet_responses/PQ4_Y.xls", 'Bx[Gauss]')
z = read_xyz("magnet_responses/PQ4_Z.xls", 'By[Gauss]')
bi = read_bi("magnet_responses/PQ4_BI.xls")


#all hardcoded for PQ4
polelen = 3000 #the length of the pole of the magnet in mm
fringelen = 300. #length of the fringe field to include outside of the magnet face in mm
nzpoints = 50 #how many z positions to include in the half magnet
zref = 1370 #z position at which the x and y scans were made
zoffset =  1782 - 3622/2.
fieldmap_current = 1100
real_current = 354



#get the field strength scaling to apply to the x, y, z fieldmaps
field_scale = interpolate_point(bi, real_current, '[A]', 'By[Gauss]')/interpolate_point(bi, fieldmap_current, '[A]', 'By[Gauss]')

x['By[Gauss]'] = field_scale*x['By[Gauss]']
y['Bx[Gauss]'] = field_scale*y['Bx[Gauss]']
z['By[Gauss]'] = field_scale*z['By[Gauss]']


x['By[Gauss]'] = x['By[Gauss]']*np.sign(x['X[mm]'])

x = include_sign(x, 'X[mm]', 'By[Gauss]', 1.)
y = include_sign(y, 'Y[mm]', 'Bx[Gauss]', 1.)


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

pd.set_option('display.max_rows', 500)
#print(zscaling)

outf = open("magnet_responses/QPQ4.dat", "w")

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
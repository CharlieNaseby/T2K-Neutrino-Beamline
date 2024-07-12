import numpy as np

g=0.0535 # T/cm

x = np.linspace(-2, 2, 9)
y = np.linspace(-2, 2, 9)
z = np.linspace(-150, 150, 51)

z = z + (1782-3622/2.)/10.

outf = open("magnet_responses/QPQ4_perfect_quad_coarse.dat", "w")

outf.write('xmin> '+str(x.min())+'\n')
outf.write('xmax> '+str(x.max())+'\n')
outf.write('nx> '+str(len(x))+'\n')

outf.write('ymin> '+str(y.min())+'\n')
outf.write('ymax> '+str(y.max())+'\n')
outf.write('ny> '+str(len(y))+'\n')

outf.write('zmin> '+str(z.min())+'\n')
outf.write('zmax> '+str(z.max())+'\n')
outf.write('nz> '+str(len(z))+'\n')

outf.write('! X \t Y \t Z \t Fx \t Fy \t Fz'+'\n')
for zpoint in z:
    for ypoint in y:
        for xpoint in x:
#            if(zpoint<z.min()+50 or zpoint>z.max()-50):
#                outf.write(str(xpoint)+' '+str(ypoint)+' '+str(zpoint)+' 0.0 0.0 0.0\n')
#            else:
            outf.write(str(xpoint)+' '+str(ypoint)+' '+str(zpoint)+' '+str(-g*ypoint)+' '+str(-g*xpoint)+' 0.0\n')
outf.close()

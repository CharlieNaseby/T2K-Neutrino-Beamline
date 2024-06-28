import numpy as np

n = 10

x = np.linspace(-10, 10, 11)
y = np.linspace(-10, 10, 11)
z = np.linspace(-200, 200, 201)

z = z + (1782-3622/2.)/10.

outf = open("magnet_responses/PQ4_dipole.dat", "w")

outf.write('xmin> -10\n')
outf.write('xmax> 10\n')
outf.write('nx> 11\n')

outf.write('ymin> -10\n')
outf.write('ymax> 10\n')
outf.write('ny> 11\n')

outf.write('zmin> '+str(z.min())+'\n')
outf.write('zmax> '+str(z.max())+'\n')
outf.write('nz> 201\n')

outf.write('! X \t Y \t Z \t Fx \t Fy \t Fz'+'\n')
for zpoint in z:
    for ypoint in y:
        for xpoint in x:
            if(zpoint<z.min()+50 or zpoint>z.max()-50):
                outf.write(str(xpoint)+' '+str(ypoint)+' '+str(zpoint)+' 0.0 0.0 0.0\n')
            else:
                outf.write(str(xpoint)+' '+str(ypoint)+' '+str(zpoint)+' 0.0 100.0 0.0\n')
outf.close()

import numpy as np
from sympy import Eq, solve
from sympy.abc import x, y
import matplotlib.pyplot as plt

def print_soln(beamwidth, emitt):
    sol = solve([ Eq(x - 2*1.684*y +1.684**2*(1+y**2)/x, beamwidth[0]**2/emitt),
                  Eq(x - 2*8.302*y +8.302**2*(1+y**2)/x, beamwidth[1]**2/emitt) ])
    print(sol)
    return sol


#beamwidth_x = 0.5*np.array([3.757336e-3, 4.242815e-3])
beamwidth_x = 0.5*np.array([3.751832e-3, 4.164317e-3])

emittx=0.084015e-6

#beamwidth_y = 0.5*np.array([1.258179e-3, 1.886464e-3])
beamwidth_y = 0.5*np.array([1.263700e-3, 1.868715e-3])
emitty=0.0695782e-6 #m rad for y


print('X solutions')
solx = print_soln(beamwidth_x, emittx)
beta = float(solx[0][x])
alpha = float(solx[0][y])

print("predicted x beam width at ssem1 (mm): ", 2000*(emittx*(beta -2*1.684*alpha + 1.684**2*(1+alpha**2)/beta))**0.5)
print("predicted x beam width at ssem2 (mm): ", 2000*(emittx*(beta -2*8.302*alpha + 8.302**2*(1+alpha**2)/beta))**0.5)


print('Y solutions')
soly = print_soln(beamwidth_y, emitty)

beta = float(soly[0][x])
alpha = float(soly[0][y])
gamma = (1+alpha**2)/beta


print("predicted y beam width at ssem1 (mm): ", 2000*(emitty*(beta -2*1.684*alpha + 1.684**2*(1+alpha**2)/beta))**0.5)
print("predicted y beam width at ssem2 (mm): ", 2000*(emitty*(beta -2*8.302*alpha + (8.302**2)*(1+alpha**2)/beta))**0.5)

s = 8.302#np.linspace(0, 10, 100)
width = beta -2*s*alpha + s**2*gamma
width = np.array(width)
width *= emitty

width = np.sqrt(width)
print(2000*width)
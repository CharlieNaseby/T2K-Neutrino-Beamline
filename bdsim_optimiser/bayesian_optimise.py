import Interface
import numpy as np



if __name__=='__main__':
    pars = [1.0, 1.0]
    target = [[0, 0, 1.0, 1.0]] #target position x, position y, width x, width y

    inter = Interface.Interface(pars, target)
    inter.SetInitialValues(pars)
    chisq = inter.fcn(pars)
    print(f'got a return value from the chisq evaluator of {chisq}')
    print(inter.fcn([1.0, 0.5]))
    #inter.bds.BeamOn(100, pars)
    #inter.GetBeamPars()
#    print(inter.fcn(testparvals))

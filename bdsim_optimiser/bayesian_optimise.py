import Interface
import numpy as np




if __name__=='__main__':
    testparvals = [i for i in range(0,21)]
    inter = Interface.Interface("ssem_data/run0910217_gen.root", "../survey/unoptimised.gmad", 21, 11, 10)
#    print(inter.fcn(testparvals))

import pandas as pd
from prody import *
from pylab import *
ion()
import numpy as np
from sklearn.preprocessing import scale

def cal_dynamics(path,name):
    ampar_ca = parsePDB(path, subset='ca')
    anm_ampar = ANM('AMPAR MT')
    anm_ampar.buildHessian(ampar_ca)
    anm_ampar.calcModes('all')
    prs_mat, effectiveness, sensitivity = calcPerturbResponse(anm_ampar)
    dfi=calcDynamicFlexibilityIndex(anm_ampar,ampar_ca,"all",norm="True")
    gnm_ampar = GNM(name)
    gnm_ampar.buildKirchhoff(ampar_ca)
    gnm_ampar.calcModes('all')
    msf=calcSqFlucts(gnm_ampar)
    stiff=calcMechStiff(anm_ampar,ampar_ca)
    newstiff=np.mean(stiff,1)
    dyn_data = np.vstack((effectiveness,
                         sensitivity,
                         msf,
                         dfi,
                         newstiff))
    return dyn_data

def dynamics_dat(name,position,WT_PDB):
    '''
    calculate dynamics using prody package, including five features:
    Effectiveness, sensitivity, stiffness, DFI and MSF
    you should prepare your protein name fot this function
    :return: a list with all dynamics calculated for the given position based on your wild type PDB
    '''
    dyn_data = cal_dynamics(WT_PDB,name)
    print("Kinetic features calculating...", end=' ')
    dyn_data = np.transpose(dyn_data)
    new_dyn = dyn_data[[x - 1 for x in position]]
    print("Done")
    return new_dyn


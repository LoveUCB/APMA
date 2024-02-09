import pandas as pd
from prody import *
from pylab import *
ion()
import numpy as np
from .. import WT_PDB
from .. import position
from .. import Protein_name

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
    dyn_data = np.vstack((effectiveness,sensitivity,msf,dfi,newstiff))
    return dyn_data

def dynamics_dat(name,position):
    dyn_data = cal_dynamics(WT_PDB,name)
    dyn_data = np.transpose(dyn_data)
    new_dyn = pd.DataFrame(columns=dyn_data.columns)
    for pos in position:
        selected_rows = dyn_data.iloc[pos - 1]
        new_dyn = pd.concat([new_dyn, selected_rows])
    return new_dyn


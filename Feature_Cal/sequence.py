from prody import *
from pylab import *
import pandas as pd
ion()
import numpy as np
import io
import sys
MSA_data = "/Users/wangjingran/Desktop/APMA/data/query_msa.fasta"
def cal_coevolution(path,position):
    msa = parseMSA(path)
    msa_refine = refineMSA(msa, label='Input_seq', rowocc=0.8, seqid=0.98)
    # showMSAOccupancy(msa_refine, occ='res')
    # coevolution
    print("Coevolution Calculating...", end = " ")
    MI = buildMutinfoMatrix(msa_refine)
    #showMutinfoMatrix(MI)
    MI = [sum(sublist)/len(sublist) for sublist in MI]
    print("Done")
    new_MI = []
    for pos in position:
        selected_rows = MI[pos - 1]
        new_MI.append(selected_rows)
    return new_MI

def cal_entropy(path,position):
    msa = parseMSA(path)
    msa_refine = refineMSA(msa, label='Input_seq', rowocc=0.8, seqid=0.98)
    #showMSAOccupancy(msa_refine, occ='res')
    # entropy
    print("Entropy Calculating...", end=" ")
    SI = calcShannonEntropy(msa_refine)
    print("Done")
    new_SI = []
    for pos in position:
        selected_rows = SI[pos - 1]
        new_SI.append(selected_rows)
    return new_SI


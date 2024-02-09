from prody import *
from pylab import *
ion()
import numpy as np

msa = parseMSA('query_msa.fasta')
msa_refine = refineMSA(msa, label='Input_seq', rowocc=0.8, seqid=0.98)
showMSAOccupancy(msa_refine, occ='res')

#conservation
showShannonEntropy(msa_refine)

# coevolution
MI = buildMutinfoMatrix(msa_refine)
showMutinfoMatrix(MI)

# entropy
Si = calcShannonEntropy(msa_refine)

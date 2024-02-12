import os
import re
import pandas as pd
print("Produced By Wang JingRan. All rights reserved.")
'''
WT_PDB = input("Please provide your route to Wild Type PDB")
Mut_PDB = input("Please provide your route ro Mutation PDB")
Protein_name = input("Please provide your protein name")
MSA_data = input("Please provide your MSA file")
description_data = input("Please provide your description file")
'''
WT_PDB = "C:/Users/33385/Desktop/APMA/APMA/data/alphafoldpten.pdb"
Mut_PDB = "C:/Users/33385/Desktop/APMA/APMA/data/PDB/Mutant"
Protein_name = "pten"
MSA_data = "C:/Users/33385/Desktop/APMA/APMA/data/query_msa.fasta"
file_path = "C:/Users/33385/Desktop/APMA/APMA/data/position.txt"

phenotype_list = []
site_list = []
with open(file_path, 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        phenotype_list.append(columns[0])
        site_list.append(columns[1])
category = phenotype_list[1:]
set_category = list(set(category))
position = site_list[1:]
position = [int(num) for num in position]


# 氨基酸网络
'''
from Feature_Cal.AAWeb import AAWEB
import os
relative_path = "data/dssp-3.0.0.exe"
absolute_path = os.path.abspath(relative_path)
absolute_path = absolute_path.replace("\\", "/")
print(absolute_path)
for i in set_category:
    AAWEB(absolute_path,i,category,Mut_PDB,WT_PDB)
'''
from Feature_Cal.AAWeb import data_AAW_gener
AAWeb_data = data_AAW_gener(position,category)

# 计算熵和保守性
from Feature_Cal.sequence import cal_entropy
from Feature_Cal.sequence import cal_coevolution
SI = cal_entropy(MSA_data,position)
MI = cal_coevolution(MSA_data,position)
# 计算弹性网络参数
from Feature_Cal.prody_cal import dynamics_dat
dynamics = dynamics_dat(Protein_name, position,WT_PDB)
df_all = pd.DataFrame()

df_all["Disease"] = category
df_all["Site"] = position

df_all["Co.evolution"] = MI
df_all["Entropy"] = SI

df_all["Betweenness"] = [sublist[0] for sublist in AAWeb_data]
df_all["Closeness"] = [sublist[1] for sublist in AAWeb_data]
df_all["Degree"] = [sublist[2] for sublist in AAWeb_data]
df_all["Eigenvector"] = [sublist[3] for sublist in AAWeb_data]
df_all["Clustering.coefficient"] = [sublist[4] for sublist in AAWeb_data]

df_all["Effectiveness"] = [sublist[0] for sublist in dynamics]
df_all["Sensitivity"] = [sublist[1] for sublist in dynamics]
df_all["MSF"] = [sublist[2] for sublist in dynamics]
df_all["DFI"] = [sublist[3] for sublist in dynamics]
df_all["Stiffness"] = [sublist[4] for sublist in dynamics]

df_all.to_csv("data/paras_old.txt", sep='\t',index=False)

print("..Machine Learning Starting...")
from ML.figure import plot_spearman
plot_spearman("data/paras.txt","Figure")
from ML import ML_Build



ML_Build()

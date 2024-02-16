import os
import re
import pandas as pd


def APMA():
    print("==================================")
    print("= Auto Protein Mutation Analyzer =")
    print("==================================")
    Protein_name = input("Please provide your protein name")
    file_path = input("Please provide your description file")
    FoldX = input("Please provide your route to FoldX")
    WT_PDB = input("Please provide your route to Wild Type PDB")
    Mut_PDB = FoldX
    MSA_data = input("Please provide your MSA file")
    '''
    FoldX = "C:/Users/33385/Desktop/FoldX"
    WT_PDB = "C:/Users/33385/Desktop/data/alphafoldpten.pdb"
    Mut_PDB = FoldX
    Protein_name = "pten"
    MSA_data = "C:/Users/33385/Desktop/data/query_msa.fasta"
    file_path = "C:/Users/33385/Desktop/data/position.txt"
    '''
    phenotype_list = []
    site_list = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(re.findall(r'\d+', columns[1])) == 1:
                site_list.append(re.findall(r'\d+', columns[1])[0])
            phenotype_list.append(columns[0])
    category = phenotype_list
    set_category = list(set(category))
    position = site_list
    position = [int(num) for num in position]

    # 使用foldx构建突变体的pdb
    from mutation.FoldX import run_FoldX
    from mutation.FoldX import get_total_energy
    run_FoldX(FoldX,WT_PDB,file_path)
    tte = get_total_energy(FoldX,WT_PDB)


    # 氨基酸网络
    from Feature_Cal.AAWeb import AAWEB
    relative_path = "data/dssp-3.0.0.exe"
    absolute_path = os.path.abspath(relative_path)
    absolute_path = absolute_path.replace("\\", "/")
    # print(absolute_path)
    print("Calculating Amino Acid Web Features", end = " ")
    for i in set_category:
        AAWEB(absolute_path,i,category,Mut_PDB,WT_PDB,"data/AAWeb")
    from Feature_Cal.AAWeb import data_AAW_gener
    AAWeb_data = data_AAW_gener(position,category)
    print("Done")
    # 计算熵和保守性
    from Feature_Cal.sequence import cal_entropy
    from Feature_Cal.sequence import cal_coevolution
    SI = cal_entropy(MSA_data,position)
    MI = cal_coevolution(MSA_data,position)
    # 计算蛋白质的相对可及表面积
    from Feature_Cal.DSSP_RASA import DSSP_RASA
    RASA = DSSP_RASA(position,"data/alphafoldpten.pdb")
    # 计算弹性网络参数
    from Feature_Cal.prody_cal import dynamics_dat
    dynamics = dynamics_dat(Protein_name, position,WT_PDB)
    df_all = pd.DataFrame()


    df_all["Disease"] = category
    df_all["Site"] = position

    df_all["Co.evolution"] = MI
    df_all["Entropy"] = SI
    df_all["RASA"] = RASA
    df_all["Total Energy"] = tte

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

if __name__ == "__main__":
    APMA()
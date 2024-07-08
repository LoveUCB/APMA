# 优化后的函数
# 调用两个线程同时计算，更加快速

import os
import re
import pandas as pd
import subprocess
import time
import threading


def APMA(WT_PDB, Protein_name, file_path, MSA_data = "/home/wangjingran/APMA/data/query_msa.fasta", FoldX = "/home/wangjingran/APMA/FoldX"):
    '''
    APMA core control panel, two threads to conduct all works.

    Running route:
    - route 1: Mutations -> FoldX -> NACEN
    - route 2: Sequence -> Blastp -> Clustal Omega -> Rate4site
    - route 1 & 2 generate 15 parameters
    - stacking model is built based on 15 features
    - model explanation based on shap

    Features:
    1. Entropy: measures the mutation frequency on the position
    2. Coevolution: coevolution on the position
    3. Conservation: conservation score based on rate4site
    4. ddG: energy change based on FoldX
    5. RASA: single residue's relative exposed area
    6. Polarity: mutation's global effect on polarity sum(MT_Polarity - WT_Polarity)
    7. Hydrophobicity: mutation's global effect on hydrophobicity sum(MT_Hydrophobicity - WT_Hydrophobicity)
    8. Betweenness: mutations' combined global effect on betweenness on a residue rowmeans(matrix(MT_betweenness) - WT_betweenness)[position]
    9. Closeness: mutations' combined global effect on closeness on a residue rowmeans(matrix(MT_closeness) - WT_closeness)[position]
    10. Eigenvector: mutations' combined global effect on eigenvector on a residue rowmeans(matrix(MT_eigenvector) - WT_eigenvector)[position]
    11. Effectiveness: dynamic network features based on ProDy
    12. Sensitivity: dynamic network features based on ProDy
    13. DFI: dynamic network features based on ProDy
    14. MSF: dynamic network features based on ProDy
    15. Sensitivity: dynamic network features based on ProDy

    Links:
    - ProDy: http://prody.csb.pitt.edu
    - FoldX: https://foldxsuite.crg.eu
    - NACEN: http://sysbio.suda.edu.cn/NACEN
    - Clustal Omega: http://www.clustal.org/omega
    - Rate4site: https://www.tau.ac.il/~itaymay/cp/rate4site.html
    - Blast: https://blast.ncbi.nlm.nih.gov
    - SHAP: https://shap.readthedocs.io/en/latest
    - Uniref50: https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50

    Parameters:
    - WT_PDB: wild type PDB file
    - Protein_name: submitted protein name
    - file_path: mutation file path
    - MSA_data: MSA data path | default: /home/wangjingran/APMA/data/query_msa.fasta
    - FoldX: FoldX folder path | default: /home/wangjingran/APMA/FoldX
    '''
 
    from Feature_Cal.Blast_MSA import extract_sequence_from_pdb
    from Feature_Cal.Blast_MSA import blast_search
    from Feature_Cal.Blast_MSA import run_clustal
    # Fetch basic informations
    # Phenotypes
    phenotype_list = []
    # Positions
    site_list = []
    # how to mutate e.g. M3V
    mutation_list = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(re.findall(r'\d+', columns[1])) == 1:
                site_list.append(re.findall(r'\d+', columns[1])[0])
            phenotype_list.append(columns[0])
            mutation_list.append(columns[1][:1] + columns[1][2:])

    category = phenotype_list
    # groups of phenotypea
    set_category = list(set(category))
    position = site_list
    position = [int(num) for num in position]
    # extract sequence from the pdb file
    pdb_sequences = extract_sequence_from_pdb(f'/home/wangjingran/APMA/data/{Protein_name}.pdb')
    protein_sequence = ''.join(pdb_sequences)
    sequence = protein_sequence

##############################################################################################################################
    def part_sequence():
        import time
        global Consurf_Scores

        # run blast
        # fetch the sequence of the pdb
        pdb_sequences = extract_sequence_from_pdb(f'/home/wangjingran/APMA/data/{Protein_name}.pdb')
        protein_sequence = ''.join(pdb_sequences)
        sequence = protein_sequence

        # Perform BLAST search and save results to a file
        # max 5 tries
        max_try_for_blast = 5
        current_try_for_blast = 0
        while current_try_for_blast < max_try_for_blast:
            current_try_for_blast += 1
            try:
                # The database is uniref50
                output_file = "/home/wangjingran/APMA/data/blast_results.fasta"
                print(f"[INFO] BLAST Search Started {current_try_for_blast} time")
                blast_search(sequence,'/home/wangjingran/prdatabase/uniref50', output_file)
                print(f"[INFO] BLAST Search success")
                break
            except Exception as e:
                print(f"[ERROR] Blast search failed {current_try_for_blast} times, {5 - current_try_for_blast} remaining")
                print(f"[ERROR] {e}")
                time.sleep(30)
        else:
            print("[ERROR] BLAST search failed after multiple tries.")
        
        # the fasta file from blastp
        with open("/home/wangjingran/APMA/data/blast_results.fasta", "r") as f:
            sequence_blast = []
            s_lines = f.readlines()
            for i in s_lines:
                if i.startswith(">") or i == "\n":
                    pass
                else:
                    sequence_blast.append(i)
            sequence_blast = list(set(sequence_blast))
            import random
            # if seq > 200, select 200 sequences randomly
            if len(sequence_blast) > 200:
                random_numbers = random.sample(range(1, len(sequence_blast)), 200)
                sequence_blast = [sequence_blast[i] for i in random_numbers]
        
        # manage the fasta file
        with open("/home/wangjingran/APMA/data/blast_results.fasta", "w") as f:
            f.write(">Input_Seq" + "\n")
            f.write(sequence + "\n")
            for i in range(len(sequence_blast)):
                f.write(">sequence" + str(i + 1) + "\n")
                f.write(sequence_blast[i])
        del sequence_blast

        # 5 max tries for msa
        max_try_for_cl = 5
        current_try_for_cl = 0
        while current_try_for_cl < max_try_for_cl:
            current_try_for_cl += 1
            try:
                # the input fasta file
                input_fasta = "/home/wangjingran/APMA/data/blast_results.fasta"
                # the output fasta file
                output_fasta = "/home/wangjingran/APMA/data/query_msa.fasta"
                # run the msa
                print(f"[INFO] ...MSA started {current_try_for_cl} time...")
                run_clustal(input_fasta, output_fasta)
                with open("/home/wangjingran/APMA/data/query_msa.fasta", 'r') as f:
                    lines = f.readlines()
                lines[0] = '>Input_seq\n'
                
                with open("/home/wangjingran/APMA/data/query_msa.fasta", 'w') as f:
                    f.writelines(lines)
                # 成功就跳出
                print("[INFO] MSA success")
                break
            except Exception as e:
                print(f"[ERROR] Clustal run failed {current_try_for_cl} times, {5 - current_try_for_cl} remaining")
                print(f"[ERROR] {e}")
                time.sleep(5)
        else:
            print("Error: MSA failed after multiple tries.")
##############################################################################################################################
        # use rate4site to score each site of the protein
        # max 5 for rate4site
        max_try_for_ra = 5
        current_try_for_ra = 0
        while current_try_for_ra < max_try_for_ra:
            current_try_for_ra += 1
            try:
                import time
                print(f"[INFO] rate4site started {current_try_for_ra} time")
                # from useless.rate4site import run_rate4site
                # time.sleep(30)
                # run_rate4site("/home/wangjingran/APMA/data/query_msa.fasta", "/home/wangjingran/APMA/data/score.txt")
                
                # excute the rate4site
                process = subprocess.Popen("rate4site -s /home/wangjingran/APMA/data/query_msa.fasta -o /home/wangjingran/APMA/data/score.txt", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                output, error = process.communicate()
                print("[INFO] Output:", output.decode().strip())
                print("[ERROR] ERROR:", error.decode().strip())

                # fetch consurf_score from rate4site
                Consurf_Score = []
                f = open("/home/wangjingran/APMA/data/score.txt","r")
                all = f.readlines()
                for i in range(len(all)):
                    if i in [0,1,2,3,4,5,6,7,8,9,10,11,12,len(all)-2,len(all)-1,len(all)]:
                        pass
                    else:
                        consurf = all[i].split()[2]
                        Consurf_Score.append(float(consurf))
                f.close()
                Consurf_Scores = []
                for i in position:
                    Consurf_Scores.append(Consurf_Score[i-1])
                print("[INFO] rate4site success")
                break
            except Exception as e:
                print(f"[ERROR] rate4site run failed {current_try_for_ra} times, {5 - current_try_for_ra} remaining")
                print(f"[ERROR] {e}")
                time.sleep(5)
        else:
            print("[ERROR] rate4site failed after multiple tries.")
        
        if not Consurf_Scores:
            raise ValueError("[ERROR] rate4site failed")



##############################################################################################################################
    # use foldx to generate pdb file and ddG
    def part_FoldX():
        global tte
        global AAWeb_data
        global dNW_data
        Mut_PDB = FoldX
        from mutation.FoldX import run_FoldX
        from mutation.FoldX import get_total_energy
        run_FoldX(FoldX,WT_PDB,file_path)
        tte = get_total_energy(FoldX,WT_PDB)


        # Amino Acid Network Featuers
        from Feature_Cal.AAWeb import AAWEB
        relative_path = "/usr/bin/mkdssp"
        absolute_path = os.path.abspath(relative_path)
        absolute_path = absolute_path.replace("\\", "/")
        print("[INFO] ...Calculating Amino Acid Network Features...", end = " ")
        
        # Cluster the amino acid contact network of each group
        # Graph Theory and statistics
        # Each mutation has an impact on the population
        # Statistics in a group
        # The cumulative value of each mutation pair site

        # Problem: sensitive to the site
        # But not sensitive to mutations

        for i in set_category:
            AAWEB(absolute_path,i,category,Mut_PDB,WT_PDB,"/home/wangjingran/APMA/data/AAWeb")
        # AAWEB(absolute_path,category,Protein_name,Mut_PDB,WT_PDB,"/home/wangjingran/APMA/data",position)
        from Feature_Cal.AAWeb import data_AAW_gener
        from Feature_Cal.AAWeb import dNW_gener
        # 获取计算出来的中心性数据
        AAWeb_data = data_AAW_gener(position,category)
        dNW_data = dNW_gener()
        print("[INFO] Animo Acid Network Features done")
    
    # generate 2 threads
    # running route 2
    thread1 = threading.Thread(target=part_sequence)
    # running route 1
    thread2 = threading.Thread(target=part_FoldX)

    # start the treads
    thread1.start()
    thread2.start()

    # wait until all threads are done
    thread1.join()
    thread2.join()
##############################################################################################################################
    # calculate entropy and coevolution
    from Feature_Cal.sequence import cal_entropy
    from Feature_Cal.sequence import cal_coevolution
    SI = cal_entropy(MSA_data,position)
    MI = cal_coevolution(MSA_data,position)
##############################################################################################################################
    # calculate rasa
    from Feature_Cal.DSSP_RASA import DSSP_RASA
    RASA = DSSP_RASA(position,WT_PDB)
##############################################################################################################################
    # calculate dynamic network features
    from Feature_Cal.prody_cal import dynamics_dat
    dynamics = dynamics_dat(Protein_name, position,WT_PDB)
##############################################################################################################################
    # The paras for ML
    df_all = pd.DataFrame()

    df_all["Disease"] = category
    df_all["Site"] = position
    df_all["Mutation"] = mutation_list

    df_all["Co.evolution"] = MI
    df_all["Entropy"] = SI
    df_all["Consurf_Score"] = Consurf_Scores
    df_all["RASA"] = RASA
    df_all["ddG"] = tte

    df_all["Betweenness"] = [sublist[0] for sublist in AAWeb_data]
    df_all["Closeness"] = [sublist[1] for sublist in AAWeb_data]
    # df_all["Degree"] = [sublist[2] for sublist in AAWeb_data]
    df_all["Eigenvector"] = [sublist[2] for sublist in AAWeb_data]

    # Polarity and Hydrophobicity
    df_all["Polarity"] = [sublist[0] for sublist in dNW_data]
    df_all["Hydrophobicity"] = [sublist[1] for sublist in dNW_data]
    
    '''
    df_all["Betweenness"] = AAWeb_data[0]
    df_all["Closeness"] = AAWeb_data[1]
    df_all["Degree"] = AAWeb_data[2]
    df_all["Eigenvector"] = AAWeb_data[3]
    df_all["Clustering.coefficient"] = AAWeb_data[4]
    '''
    
    df_all["Effectiveness"] = [sublist[0] for sublist in dynamics]
    df_all["Sensitivity"] = [sublist[1] for sublist in dynamics]
    df_all["MSF"] = [sublist[2] for sublist in dynamics]
    df_all["DFI"] = [sublist[3] for sublist in dynamics]
    df_all["Stiffness"] = [sublist[4] for sublist in dynamics]

    # save into txt file
    df_all.to_csv("/home/wangjingran/APMA/data/paras.txt", sep='\t',index=False)
    df_all.to_csv("/home/wangjingran/APMA/Outcome/paras.txt",sep = '\t', index=False)
############################################################################################################################## 
    # ML module for the features
    print("[INFO] ...Machine Learning Starting...")
    from ML import ML_Build
    ML_Build(category)

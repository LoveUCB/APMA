# 优化后的函数
# 调用两个线程同时计算，更加快速

import os
import re
import pandas as pd
import subprocess
import time
import threading



def APMA(WT_PDB, Protein_name, file_path,MSA_data = "/home/wangjingran/APMA/data/query_msa.fasta", FoldX = "/home/wangjingran/APMA/FoldX"):
    '''
    APMA 机器学习辅助的整合蛋白质突变模型
    作者：王景然
    单位：苏州大学医学院生物信息
    参数：
        WT_PDB: 野生型蛋白质文件的位置
        Protein_name: 蛋白质的名称,默认为蛋白质结构文件的文件名
        file_path: 记录突变的文件位置
    '''
 
    from Feature_Cal.Blast_MSA import extract_sequence_from_pdb
    from Feature_Cal.Blast_MSA import blast_search
    from Feature_Cal.Blast_MSA import run_clustal
    # 基本信息
    # 标签:属于哪一组
    phenotype_list = []
    # 突变发生在蛋白质的哪个位点
    site_list = []
    # 突变的形式 如M3V
    mutation_list = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(re.findall(r'\d+', columns[1])) == 1:
                site_list.append(re.findall(r'\d+', columns[1])[0])
            phenotype_list.append(columns[0])
            mutation_list.append(columns[1][:1] + columns[1][2:])

    category = phenotype_list
    # 获取数据有几组
    set_category = list(set(category))
    position = site_list
    position = [int(num) for num in position]
    # 获取蛋白质的全序列
    pdb_sequences = extract_sequence_from_pdb(f'/home/wangjingran/APMA/data/{Protein_name}.pdb')
    protein_sequence = ''.join(pdb_sequences)
    sequence = protein_sequence

    # 先进行blast
    # 获取蛋白质的全序列
    pdb_sequences = extract_sequence_from_pdb(f'/home/wangjingran/APMA/data/{Protein_name}.pdb')
    protein_sequence = ''.join(pdb_sequences)
    sequence = protein_sequence

    # Perform BLAST search and save results to a file
    # 搜索可能会失败，设置最多五次
    max_try_for_blast = 3
    current_try_for_blast = 0
    while current_try_for_blast < max_try_for_blast:
        current_try_for_blast += 1
        try:
            output_file = "/home/wangjingran/APMA/data/blast_results.fasta"
            print(f"BLAST Search Started {current_try_for_blast} time")
            blast_search(sequence, output_file)
            print(f"BLAST Search success")
            break
        except Exception as e:
            print(f"Blast search failed {current_try_for_blast} times, {5 - current_try_for_blast} remaining")
            print(f"Error: {e}")
            time.sleep(30)
    else:
        print("Error: BLAST search failed after multiple tries.")
    
    # 输入的FASTA文件，这里假设你已经有了一些同源序列的FASTA文件
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
        # 这里加上一个判断，如果少于了200个就把所有的都选上去
        if len(sequence_blast) > 200:
            random_numbers = random.sample(range(1, len(sequence_blast)), 200)
            sequence_blast = [sequence_blast[i] for i in random_numbers]
    
    with open("/home/wangjingran/APMA/data/blast_results.fasta", "w") as f:
        f.write(">Input_Seq" + "\n")
        f.write(sequence + "\n")
        for i in range(len(sequence_blast)):
            f.write(">sequence" + str(i + 1) + "\n")
            f.write(sequence_blast[i])
    del sequence_blast
##############################################################################################################################
    def part_sequence():
        import time
        global Consurf_Scores
        max_try_for_cl = 5
        current_try_for_cl = 0
        while current_try_for_cl < max_try_for_cl:
            current_try_for_cl += 1
            try:
                input_fasta = "/home/wangjingran/APMA/data/blast_results.fasta"
                # 输出的FASTA文件，用于保存比对结果
                output_fasta = "/home/wangjingran/APMA/data/query_msa.fasta"
                
                # 运行多序列比对
                print(f"MSA started {current_try_for_cl} time")
                run_clustal(input_fasta, output_fasta)
                with open("/home/wangjingran/APMA/data/query_msa.fasta", 'r') as f:
                    lines = f.readlines()
                lines[0] = '>Input_seq\n'
                
                with open("/home/wangjingran/APMA/data/query_msa.fasta", 'w') as f:
                    f.writelines(lines)
                # 成功就跳出
                print("MSA success")
                break
            except Exception as e:
                print(f"Clustal run failed {current_try_for_cl} times, {5 - current_try_for_cl} remaining")
                print(f"Error: {e}")
                time.sleep(5)
        else:
            print("Error: MSA failed after multiple tries.")
##############################################################################################################################
        # use rate4site to score each site of the protein
        max_try_for_ra = 5
        current_try_for_ra = 0
        while current_try_for_ra < max_try_for_ra:
            current_try_for_ra += 1
            try:
                import time
                print(f"rate4site started {current_try_for_ra} time")
                # from useless.rate4site import run_rate4site
                # time.sleep(30)
                # run_rate4site("/home/wangjingran/APMA/data/query_msa.fasta", "/home/wangjingran/APMA/data/score.txt")
                
                process = subprocess.Popen("rate4site -s /home/wangjingran/APMA/data/query_msa.fasta -o /home/wangjingran/APMA/data/score.txt", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                output, error = process.communicate()
                print("Output:", output.decode().strip())
                print("Error:", error.decode().strip())

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
                print("rate4site success")
                break
            except Exception as e:
                print(f"rate4site run failed {current_try_for_ra} times, {5 - current_try_for_ra} remaining")
                print(f"Error: {e}")
                time.sleep(5)
        else:
            print("Error: rate4site failed after multiple tries.")
        
        if not Consurf_Scores:
            raise ValueError("rate4site failed")
##############################################################################################################################
    # 使用foldx构建突变体的pdb
    def part_FoldX():
        global tte
        global AAWeb_data
        Mut_PDB = FoldX
        from mutation.FoldX import run_FoldX
        from mutation.FoldX import get_total_energy
        run_FoldX(FoldX,WT_PDB,file_path)
        tte = get_total_energy(FoldX,WT_PDB)


        # 氨基酸网络
        from Feature_Cal.AAWeb import AAWEB
        relative_path = "/usr/bin/mkdssp"
        absolute_path = os.path.abspath(relative_path)
        absolute_path = absolute_path.replace("\\", "/")
        print("Calculating Amino Acid Web Features...", end = " ")
        
        # 对每一组的氨基酸接触网络分组聚类
        # 图论与统计
        # 每一个突变对总体均会造成影响
        # 统计在一个组中
        # 每个突变对位点的累加值

        # 有问题：对位点敏感
        # 但是对突变不敏感
        for i in set_category:
            AAWEB(absolute_path,i,category,Mut_PDB,WT_PDB,"/home/wangjingran/APMA/data/AAWeb")
        # AAWEB(absolute_path,category,Protein_name,Mut_PDB,WT_PDB,"/home/wangjingran/APMA/data",position)
        from Feature_Cal.AAWeb import data_AAW_gener
        # 获取计算出来的中心性数据
        AAWeb_data = data_AAW_gener(position,category)
        print("Done")
    
    # 创建线程
    # 将两个最耗时间的分开计算
    thread1 = threading.Thread(target=part_sequence)
    thread2 = threading.Thread(target=part_FoldX)

    # 启动线程
    thread1.start()
    thread2.start()

    # 等待两个线程都执行完毕
    thread1.join()
    thread2.join()
##############################################################################################################################
    # 计算熵和保守性
    from Feature_Cal.sequence import cal_entropy
    from Feature_Cal.sequence import cal_coevolution
    SI = cal_entropy(MSA_data,position)
    MI = cal_coevolution(MSA_data,position)
##############################################################################################################################
    # 计算蛋白质的相对可及表面积
    from Feature_Cal.DSSP_RASA import DSSP_RASA
    RASA = DSSP_RASA(position,WT_PDB)
##############################################################################################################################
    # 计算弹性网络参数
    from Feature_Cal.prody_cal import dynamics_dat
    dynamics = dynamics_dat(Protein_name, position,WT_PDB)
##############################################################################################################################
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
    df_all["Degree"] = [sublist[2] for sublist in AAWeb_data]
    df_all["Eigenvector"] = [sublist[3] for sublist in AAWeb_data]
    df_all["Clustering.coefficient"] = [sublist[4] for sublist in AAWeb_data]
    
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

# 将结果保存到paras.txt文件中
    df_all.to_csv("/home/wangjingran/APMA/data/paras.txt", sep='\t',index=False)
    df_all.to_csv("/home/wangjingran/APMA/Outcome/paras.txt",sep = '\t', index=False)
############################################################################################################################## 
    print("..Machine Learning Starting...")
    from ML import ML_Build
    ML_Build(category)

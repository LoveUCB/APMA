from APMA import APMA
try:
    APMA(
        Protein_name = "pten",
        file_path = "/home/wangjingran/APMA/data/position.txt",
        WT_PDB = "/home/wangjingran/APMA/data/alphafoldpten.pdb"
        )


    from ML.figure import plot_roc_for_disease_pairs
    plot_roc_for_disease_pairs("data/paras.txt","/home/wangjingran/APMA/Outcome/Figure/ROC/Feature")

    from ML.figure import plot_box
    plot_box("data/paras.txt","/home/wangjingran/APMA/Outcome/Figure/Box_Violin")

    from ML.figure import plot_spearman
    plot_spearman("data/paras.txt","/home/wangjingran/APMA/Outcome/Figure")

    from Email.zip import zip_folder
    zip_folder('/home/wangjingran/APMA/Outcome','/home/wangjingran/APMA/Email/APMA_outcome.zip')
    from Email.send import send_email
    send_email(['3338561620@qq.com'])
except:
    from Email.send import send_error_email
    send_error_email(['3338561620@qq.com'])


import os

def delete_files_in_directory(directory):
    # 遍历目录中的所有文件和子目录
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        # 如果是文件，则删除
        if os.path.isfile(item_path):
            os.remove(item_path)
        # 如果是目录，则递归调用该函数
        elif os.path.isdir(item_path):
            delete_files_in_directory(item_path)

delete_files_in_directory("/home/wangjingran/APMA/Outcome")
delete_files_in_directory("/home/wangjingran/APMA/data")

folder_path = '/home/wangjingran/APMA/FoldX'
files = os.listdir(folder_path)

for file_name in files:
    if file_name != 'foldx4' and file_name != 'rotabase.txt':
        file_path = os.path.join(folder_path, file_name)

        os.remove(file_path)


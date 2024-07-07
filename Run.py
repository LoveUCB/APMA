# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/APMA


"""

#############################################
### Introduction of Run Module
#
# @ This module is to control the main progress
#
#############################################


from APMA import APMA
# from . import __yourdownloadroute__
# APMA_path = __yourdownloadroute__
import os
import glob
import traceback
import datetime
import re
import requests

# fetch start time
current_datetime = datetime.datetime.now()

# print start time
print("***** APMA Start at: ", current_datetime)


def find_consecutive_numbers(input_string):
    """
    Finds all consecutive numbers in a given string and returns them as a list of strings.

    Parameters:
    input_string (str): The input string in which to find consecutive numbers.

    Returns:
    list: A list of strings representing the consecutive numbers found in the input string.
    """
    pattern = r'\d+'
    matches = re.findall(pattern, input_string)
    return matches


def download_alphafold_structure(uniprot_id, output_path):
    """
    Downloads the protein structure file from the AlphaFold database based on the given Uniprot ID and saves it to the specified output path.

    Parameters:
    uniprot_id (str): The Uniprot ID of the protein for which the structure file is to be downloaded.
    output_path (str): The path (including filename) where the downloaded structure file should be saved.

    Returns:
    None

    Example:
    download_alphafold_structure("P12345", "P12345.pdb")
    """
    print("[INFO] ...Downloading Structure File from AlphaFold database...")
    base_url = "https://alphafold.ebi.ac.uk/files/AF-"
    file_url = f"{base_url}{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(file_url)

    if response.status_code == 200:
        with open(output_path, 'wb') as file:
            file.write(response.content)
        print(f"[INFO] Structure file for {uniprot_id} downloaded successfully.")
    else:
        raise ValueError(f"[INFO] Failed to download structure file for {uniprot_id}. Status code: {response.status_code}")


def delete_files_in_directory(directory):
    """
    Deletes all files and directories within the specified directory.
    
    Parameters:
    directory (str): The path to the directory from which to delete files and directories.
    
    Returns:
    None
    
    This function traverses the specified directory and deletes all files.
    If a subdirectory is encountered, it recursively deletes its contents as well.
    
    Example:
    delete_files_in_directory("/home/user/directory_to_clean")
    """
    # Traverse all files and subdirectories in the specified directory
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        
        # If the item is a file, delete it
        if os.path.isfile(item_path):
            os.remove(item_path)
        
        # If the item is a directory, recursively call this function
        elif os.path.isdir(item_path):
            delete_files_in_directory(item_path)


# Define Amino Acid list
Amino_acids_list = [
    'A',  # Alanine
    'R',  # Arginine
    'N',  # Asparagine
    'D',  # Aspartic acid
    'C',  # Cysteine
    'Q',  # Glutamine
    'E',  # Glutamic acid
    'G',  # Glycine
    'H',  # Histidine
    'I',  # Isoleucine
    'L',  # Leucine
    'K',  # Lysine
    'M',  # Methionine
    'F',  # Phenylalanine
    'P',  # Proline
    'S',  # Serine
    'T',  # Threonine
    'W',  # Tryptophan
    'Y',  # Tyrosine
    'V'   # Valine
]


try:
    # check whether need to download structure from AlphaFold database
    uni_file_path = "/home/wangjingran/APMA/data/uniprot_ID.txt"
    if os.path.exists(uni_file_path):
        submit = False
    else:
        submit = True

    # Get the Uniprot ID from the user
    if submit:
        pass
    else:
        f = open('/home/wangjingran/APMA/data/uniprot_ID.txt')
        uniprot_id = f.readlines()
        count_id = 0
        for i in uniprot_id:
            if i == "\n":
                pass
            else:
                count_id += 1
        if count_id == 1:
            uniprot_id = uniprot_id[0].strip('\n')
        else:
            raise ValueError("[ERROR] We do not support multiple uniprot IDs")
        f.close()
        # Download the structure file
        download_alphafold_structure(uniprot_id, f"/home/wangjingran/APMA/data/{uniprot_id}.pdb")

    # pocess position file
    dict_res_type = {}
    f = open("/home/wangjingran/APMA/data/position.txt","r")
    all = f.readlines()
    all_new = []
    # print(all)
    for i in all:
        if i == '\n':
            pass
        else:
            string_split = i.split("\t")
            string_a = string_split[0]
            string_b = string_split[1]

            # check amino acid
            ori_amino_acid = string_b[0]
            mut_amino_acid = string_b[-2]
            if (ori_amino_acid not in Amino_acids_list) or (mut_amino_acid not in Amino_acids_list):
                raise ValueError("[ERROR] Unsupported amino acid abbreviation")
            
            # change to FoldX forat
            string_b = string_b[:1] + "A" + string_b[1:]

            # add into dict
            dict_res_type[int(find_consecutive_numbers(string_b)[0])] = string_b[0]
            FoldX_type_string = string_a + "\t" + string_b
            all_new.append(FoldX_type_string)
    f.close()


    f = open("/home/wangjingran/APMA/data/position.txt","w")
    for i in all_new:
        f.write(i)
    f.close()

    del all
    del all_new

    # fetch email
    email_list = []
    f = open("/home/wangjingran/APMA/data/email.txt")
    lines = f.readlines()
    for i in lines:
        line = i.strip("\n")
        email_list.append(line)
    f.close()

    # fetch pdb file
    def print_pdb_files(folder_path):
        # 使用 glob 模块列出文件夹中所有的 .pdb 文件
        pdb_files = glob.glob(os.path.join(folder_path, '*.pdb'))
        for i in pdb_files:
            user_pdb_file = i
        return user_pdb_file

    folder_path = '/home/wangjingran/APMA/data'
    user_pdb = print_pdb_files(folder_path)
    user_protein_name = user_pdb.split("/")[-1].rstrip(".pdb")
    # print(dict_res_type)

    # check position file
    from Feature_Cal.Blast_MSA import extract_sequence_from_pdb
    pdb_seq = extract_sequence_from_pdb(user_pdb)
    pdb_seq = ''.join(pdb_seq)
    for (position, residue) in dict_res_type.items():
        if pdb_seq[position - 1] != residue:
            raise ValueError(f"[ERROR] Position {position} in the position file does not match the PDB sequence.")
    print("[INFO] Pass checking")
    
    # check email
    for i in email_list:
        if '@' not in i:
            raise ValueError("[ERROR] Invalid email address.")

    # After passing several checks, run the main function
    #####################################################
    APMA(
    Protein_name = user_protein_name,
    file_path = "/home/wangjingran/APMA/data/position.txt",
    WT_PDB = user_pdb
    )
    #####################################################

    from ML.figure import plot_roc_for_disease_pairs
    plot_roc_for_disease_pairs("/home/wangjingran/APMA/data/paras.txt","/home/wangjingran/APMA/Outcome/Figure/ROC/Feature")

    from ML.figure import plot_box
    plot_box("/home/wangjingran/APMA/data/paras.txt","/home/wangjingran/APMA/Outcome/Figure")

    from ML.figure import plot_spearman
    plot_spearman("/home/wangjingran/APMA/data/paras.txt","/home/wangjingran/APMA/Outcome/Figure")

    from ML.figure import plot_dynamic_network
    plot_dynamic_network('/home/wangjingran/APMA/data/all_dyn_data.txt', '/home/wangjingran/APMA/data/paras.txt', '/home/wangjingran/APMA/Outcome/Figure/dynamic.pdf')

    from Email.zip import zip_folder
    zip_folder('/home/wangjingran/APMA/Outcome','/home/wangjingran/APMA/Email/APMA_outcome.zip')
    
    from Email.send import send_email
    send_email(email_list)

except Exception as e:
    print(str(e))
    traceback_info = traceback.format_exc()
    print(traceback_info)
    current_datetime = datetime.datetime.now()
    print("[INFO] APMA ends at: ", current_datetime)
    from Email.send import send_error_email
    send_error_email(email_list)


delete_files_in_directory("/home/wangjingran/APMA/Outcome")
delete_files_in_directory("/home/wangjingran/APMA/data")

folder_path = '/home/wangjingran/APMA/FoldX'
files = os.listdir(folder_path)

for file_name in files:
    
    if file_name != 'foldx4' and file_name != 'rotabase.txt' and file_name != 'foldx5' and file_name != 'molecules':
        file_path = os.path.join(folder_path, file_name)

        os.remove(file_path)


# print end time
print("[INFO] APMA ends at: ", current_datetime)

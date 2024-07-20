# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/deePheMut

"""

#############################################
### Introduction of zip module
#
# @ This module is to zip Outcome folder
#
#############################################



import zipfile
import os

import os
import zipfile

def zip_folder(folder_path, zip_name):
    """
    Compresses the contents of a folder into a ZIP file.

    Parameters:
    folder_path (str): The path to the folder to be compressed.
    zip_name (str): The name of the output ZIP file.

    Raises:
    FileNotFoundError: If the specified folder does not exist.
    """
    
    # Ensure the folder path exists
    if not os.path.exists(folder_path):
        raise FileNotFoundError("Folder does not exist")
    
    # Create a ZipFile object
    with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # Iterate over all the files and subfolders in the folder
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                # Construct the full path of the file
                file_path = os.path.join(root, file)
                # Compute the relative path of the file in the ZIP file
                relative_path = os.path.relpath(file_path, folder_path)
                # Add the file to the ZIP file
                zipf.write(file_path, relative_path)


if __name__ == "__main__":
    folder_to_zip = '/home/wangjingran/APMA/data'
    zip_file_name = '/home/wangjingran/APMA/Email/Outcome.zip'
    zip_folder(folder_to_zip, zip_file_name)
    print("success")

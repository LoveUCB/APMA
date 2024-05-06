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

    if file_name != 'foldx4' and file_name != 'rotabase.txt' and file_name != 'foldx5' and file_name != 'molecules':
        file_path = os.path.join(folder_path, file_name)

        os.remove(file_path)

print("Cleared")

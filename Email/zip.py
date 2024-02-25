import zipfile
import os

def zip_folder(folder_path, zip_name):
    # 确保文件夹路径存在
    if not os.path.exists(folder_path):
        raise FileNotFoundError("Folder does not exist")

    # 创建一个ZipFile对象
    with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # 遍历文件夹中的所有文件和子文件夹
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                # 构建文件的完整路径
                file_path = os.path.join(root, file)
                # 计算文件在压缩文件中的相对路径
                relative_path = os.path.relpath(file_path, folder_path)
                # 将文件添加到压缩文件中
                zipf.write(file_path, relative_path)

# 要压缩的文件夹路径和压缩文件的名称
if __name__ == "__main__":
    folder_to_zip = '/home/wangjingran/APMA/data'
    zip_file_name = '/home/wangjingran/APMA/Email/APMA_outcome.zip'

    # 调用函数压缩文件夹
    zip_folder(folder_to_zip, zip_file_name)
    print("success")

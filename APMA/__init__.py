import os
import re
print("Produced By Wang JingRan. All rights reserved.")
WT_PDB = input("Please provide your route to Wild Type PDB")
Mut_PDB = input("Please provide your route ro Mutation PDB")
Protein_name = input("Please provide your protein name")
def get_position(path):
    file_names = os.listdir(path)
    pattern = re.compile(r'\d+')
    position = []
    for file_name in file_names:
        # 使用正则表达式查找数字部分
        match = pattern.search(file_name)
        if match:
            position.append(match.group())
    return position
position = get_position(Mut_PDB)

def get_category(path):
    file_names = os.listdir(path)
    pattern = re.compile(r'[a-zA-Z]+')  # 匹配字母部分的正则表达式
    category = []
    for file_name in file_names:
        # 使用正则表达式查找字母部分
        match = pattern.search(file_name)
        if match:
            category.append(match.group())
    return category
category = get_category(Mut_PDB)

from Feature_Cal.prody_cal import dynamics_dat
def predict():
    dyn_data = dynamics_dat(Protein_name,position)
    dyn_data["catagory"] = category




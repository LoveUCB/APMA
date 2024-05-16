import pymol
from pymol import cmd

def visualize_protein(selection_residues, pdb_path):
    # 启动PyMOL无GUI模式
    pymol.finish_launching(['pymol', '-cq'])

    # 加载蛋白质文件
    cmd.load(pdb_path)

    first_residue_selection = f"resi {selection_residues[0]}"

    # 接着选择中间氨基酸的阿尔法碳原子
    middle_residues_selection = ' or '.join([f"resi {res} and name CA" for res in selection_residues[1:-1]])

    # 最后选择最后一个氨基酸的所有原子
    last_residue_selection = f"resi {selection_residues[-1]}"

    # 组合所有选择条件
    final_selection = f"{first_residue_selection} or {middle_residues_selection} or {last_residue_selection}"
    node_selection = 'resi ' + '+'.join(map(str, selection_residues))

    # 在PyMOL中执行选择和显示命令
    cmd.select("selected_residues", final_selection)
    cmd.show("spheres", "selected_residues")
    cmd.set("sphere_scale", 1, first_residue_selection)
    cmd.set("sphere_scale", 1, last_residue_selection)
    cmd.set("sphere_scale", 0.7, middle_residues_selection)

    # 将背景蛋白质显示为cartoon并虚化
    cmd.show('cartoon', 'all')
    cmd.set('cartoon_transparency', 0.5, 'all')

    # 调整背景蛋白质的颜色，使其看起来更淡
    cmd.color('gray90', 'all')


    # 路径展示
    # 创建一个新对象
    cmd.create("my_path", "sele")

    # 选择要连接的原子，例如：
    for j in range(len(selection_residues)):
        cmd.select(f"atom{j}", f"resi {selection_residues[j]} and name CA")


    # 在选定的原子之间创建连接
    for j in range(len(selection_residues) - 1):
        cmd.bond(f"atom{j}", f"atom{j + 1}")
        # cmd.set_bond("stick_color", "selenium", selection1=f"my_path and atom{j}", selection2=f"my_path and atom{j}")
    cmd.color("selenium", "my_path")
        

    # 显示连接路径
    cmd.show("sticks", "my_path")
    cmd.color("selenium", "my_path")

    # 为每个选中的氨基酸生成颜色
    colors = ['warmpink', 'gray50', 'skyblue']
    for i, resi in enumerate(selection_residues):
        if i == 0:
            cmd.color(colors[0], f'resi {resi}')
        elif i == len(selection_residues) - 1:
            cmd.color(colors[2], f'resi {resi}')
        else:
            pass

    # 设置视角和光线追踪参数
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 0)

    # 使用光线追踪渲染并保存图片
    print("...saving figures...")
    cmd.ray(12000, 9000)
    cmd.png('output_image.png')

    # 退出PyMOL
    cmd.quit()
    print("Success")

# 例子调用函数
if __name__ == "__main__":
    visualize_protein([1, 13, 46, 77, 89, 105, 138, 189, 204], "/home/wangjingran/APMA/md-dyn-net/pten.pdb")

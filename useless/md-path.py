import networkx as nx
import matplotlib.pyplot as plt
import re
import os
import os
from collections import defaultdict
from collections import Counter
import subprocess
from pymol_show_pathway import visualize_protein

def run_md_task(path_to_pdb, step, md_file):
    f = open("prepare_md-task.sh", "w")
    f.write("conda activate md-task\n")
    f.write("export PATH=/home/wangjingran/MD-TASK:$PATH-task\n")
    f.write(f"calc_network.py --topology {path_to_pdb} --threshold 7.0 --step {step} --generate-plots --calc-L  --lazy-load {md_file}\n")
    f.close()
    print("...Start running md-task...")
    process = subprocess.Popen(f"bash prepare_md-task.sh", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    print("Output:", output.decode().strip())
    print("Error:", error.decode().strip())


def edge_frequency(networks):
    print("...Combine networks...")
    # 创建一个字典用于存储边的出现频率
    edge_weights = defaultdict(int)
    # 遍历每个网络
    for network in networks:
        # 遍历网络中的每条边
        for edge in network.edges():
            # 将边转换为元组，并作为键来更新出现频率
            edge_weights[tuple(sorted(edge))] += 1
    return edge_weights

def combine_network(folder_path, record=False):
    dyn_frames = []
    # 读取全部的gml文件
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gml"):
            file_path = os.path.join(folder_path, file_name)
            G = nx.read_gml(file_path)
            dyn_frames.append(G)
    len_network = len(dyn_frames)
    integrated_network = nx.compose_all(dyn_frames)

    if record:
        # 计算边的出现频率
        edge_weights = edge_frequency(dyn_frames)
        # 保存边的权重到文件
        output_file = record
        with open(output_file, 'w') as f:
            for edge, weight in edge_weights.items():
                f.write(f"{edge[0]}\t{edge[1]}\t{weight / len_network}\n")

    return integrated_network






def graph_short_path(file, output, start, end, cutoff, record = False, plot = True):
    '''
    一个寻找网络两个节点之间最短路径并绘图的函数
    parameters:
        file: 记录网络节点和链接的文件
        output: 保存绘图和路径搜索的位置
    '''
    G = nx.Graph()
    edge_prob_dict = {}
    print("...Building Graph...")
    f = open(file, "r")
    all = f.readlines()
    for i in all:
        m = i.split("\t")
        a = m[0].strip("\n")
        b = m[1].strip("\n")
        l = m[2].strip("\n")
        a = re.findall(r'\d+', a)[0]
        b = re.findall(r'\d+', b)[0]
        a = str(int(a) + 1)
        b = str(int(b) + 1)
        if float(l) < cutoff:
            pass
        else:
            if int(a) < int(b):
                G.add_edge(a, b)
                edge_prob_dict[f"{a}, {b}"] = float(l)
            else:
                G.add_edge(b, a)
                edge_prob_dict[f"{b}, {a}"] = float(l)
    f.close()


    print("...Searching full shortest route in unweighted graph...")
    shortest_path = nx.all_shortest_paths(G, source=start, target=end)
    shortest_path = list(shortest_path)

    prob_list = []
    for sub_path in shortest_path:
        prob = 1
        for unit in range(len(sub_path) - 1):
            m = sub_path[unit]
            n = sub_path[unit + 1]
            if int(m) < int(n):
                pass
            else:
                m, n = n, m
            prob_p = edge_prob_dict[f"{m}, {n}"]
            prob = prob * prob_p
        prob_list.append(prob)
    max_prob = max(prob_list)
    
    shortest_list_finall = []
    for i in range(len(prob_list)):
        if prob_list[i] == max_prob:
            shortest_list_finall.append(shortest_path[i])

    n = 1
    for shortest_list_final in shortest_list_finall:
    
        # 当record是True对最短路径进行记录
        if record != False:
            f = open(f"{output}/record_route.txt", "a")
            f.write(f"from {start} to {end}: \t")
            f.write(" -> ".join(shortest_list_final) + "\n")
            f.close()
            # print(f"from position {start} to position {end}")
            print("shortest route {n}:", " -> ".join(shortest_list_final))
        else:
            # print(f"from position {start} to position {end}")
            print(f"shortest route {n}:", " -> ".join(shortest_list_final))
            

        
        # 绘图
        if plot == True:
            print("...Saving Figure...")
            pos = nx.spring_layout(G, k=0.15, seed=4572321)
            plt.figure(figsize=(20, 20))
            nx.draw_networkx_nodes(G, pos, node_size=30, node_color="#82B0D2", label=True, alpha=1)
            nx.draw_networkx_edges(G, pos, width=0.2, edge_color="gainsboro", alpha=1)
            path_edges = [(shortest_list_final[i], shortest_list_final[i + 1]) for i in range(len(shortest_list_final) - 1)]
            path_nodes = shortest_list_final
            node_colors = ['#ec4347' if node in [start, end] else 'orange' for node in path_nodes]
            node_size = [400 if node in [start, end] else 300 for node in path_nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=path_nodes, node_color=node_colors, node_size=node_size)
            shortest_path_labels = {node: node for node in path_nodes}
            nx.draw_networkx_labels(G, pos, labels=shortest_path_labels, font_size=7)
            nx.draw_networkx_edges(G, pos, edgelist=path_edges, width=1.8, edge_color='orange', arrows=True, arrowstyle='->')
            plt.axis('off')
            plt.savefig(f"{output}/path from {start} to {end}_{n}.pdf")
            plt.close()
            print(f"Figure has been saved to {output}path from {start} to {end}_{n}.pdf")
        else:
            pass
        n += 1
    return shortest_list_finall


def delete_files_with_extensions(folder_path, extensions):
    # 遍历文件夹中的所有文件
    print("...Cleaning md-task files...")
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        # 检查文件是否是一个文件而不是文件夹
        if os.path.isfile(file_path):
            # 检查文件的后缀名是否在指定的扩展列表中
            if any(file_name.endswith(ext) for ext in extensions):
                # 删除文件
                os.remove(file_path)






if __name__ == "__main__":
    # pdb_file_path = input("Your route to PDB file: ")
    # md_file = input("Your route to md file: ")
    # Step = input("Step: ")
    # start_AA = input("Start: ")
    # end_AA = input("End: ")

    pdb_file_path = '/home/wangjingran/MD-TASK/R88Q.pdb'
    md_file = '/home/wangjingran/MD-TASK/net.xtc'
    Step = 10
    start_AA = '88'
    end_AA = '915'
    edge_cutoff = 0.5
    run_md_task(pdb_file_path, int(Step), md_file)

    combined_network = combine_network('./',
                                       record="./Combined_Dyn_Net.txt")
    
    sp = graph_short_path(
                './Combined_Dyn_Net.txt', 
                './', 
                start_AA, end_AA,
                cutoff = edge_cutoff,
                plot = True
                    )
    
    # 指定文件夹路径和要删除的文件后缀名列表
    folder_path = "./"
    extensions_to_delete = [".dat", ".gml", ".graphml", ".png"]
    # 调用函数删除指定后缀名的文件
    delete_files_with_extensions(folder_path, extensions_to_delete)
    n = 1
    for i in sp:
        print(f"...PyMOL Running for route {n}...")
        visualize_protein(i, pdb_file_path, n)
        n += 1
        
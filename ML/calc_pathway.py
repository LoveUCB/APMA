import networkx as nx
import matplotlib.pyplot as plt
import re
import os
from tqdm import tqdm


def graph_short_path(file, output, start, end, record = False, plot = True, length = False, full_route = True):
    '''
    一个寻找网络两个节点之间最短路径并绘图的函数
    parameters:
        file: 记录网络节点和链接的文件
        output: 保存绘图和路径搜索的位置
    '''
    if length:
        if full_route:
            print("** full shortest route search in weighted graph **")
        else:
            print("** restricted shortest route search in weighted graph **")
    else:
        if full_route:
            print("** full shortest route search in unweighted graph **")
        else:
            print("** restricted shortest route search in unweighted graph **")
    # 初始化nx模块
    G = nx.Graph()
    if not length:
    # 打开记录网络的文件
        print("...Building Graph...", end=' ')
        f = open(file,"r")
        all = f.readlines()
        for i in all:
            m = i.split("\t")
            a = m[0].strip("\n")
            b = m[1].strip("\n")
            a = re.findall(r'\d+', a)[0]
            b = re.findall(r'\d+', b)[0]
            a = str(int(a) + 1)
            b = str(int(b) + 1)
            G.add_edge(a, b)
        f.close()
        print("Success")

        if full_route:
            print("...Searching full shortest route in unweighted graph...")
            shortest_path = nx.all_shortest_paths(G, source=start, target=end)
            shortest_path = list(shortest_path)
        else:
            print("...Searching restricted shortest route in unweighted graph...")
            shortest_path = nx.shortest_path(G, source=start, target=end)
            shortest_path = list(shortest_path)

    else:
        print("...Building Graph...", end=' ')
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
            l = re.findall(r'\d+', l)[0]
            G.add_edge(a, b, weight = int(l))
        f.close()
        print("Success")

        if full_route:
            print("...Searching full shortest route in weighted graph...")
            shortest_path = nx.all_shortest_paths(G, source=start, target=end, weight='weight')
            shortest_path = list(shortest_path)
        else:
            print("...Searching restricted shortest route in weighted graph...")
            shortest_path = nx.shortest_path(G, source=start, target=end, weight='weight')
            shortest_path = list(shortest_path)
    

    # 当record是True对最短路径进行记录
    if full_route:
        if record != False:
            f = open(f"{output}/record_route.txt", "a")
            for i in shortest_path:
                f.write(f"from {start} to {end}: \t")
                f.write(" -> ".join(i) + "\n")
            f.close()
            print(f"from position {start} to position {end}")
            for i in shortest_path:
                print("shortest route:", " -> ".join(i))
        else:
            print(f"from position {start} to position {end}")
            print(f"There are {len(shortest_path)} shortest route detected")
            count_route = 1
            for i in shortest_path:
                print(f"shortest route num {count_route}:", " -> ".join(i))
                count_route += 1
    
    else:
        if record != False:
            f = open(f"{output}/record_route.txt", "a")
            f.write(f"from {start} to {end}: \t")
            f.write(" -> ".join(shortest_path) + "\n")
            f.close()
            print(f"from position {start} to position {end}")
            print("shortest route:", " -> ".join(shortest_path))
        else:
            print(f"from position {start} to position {end}")
            print(f"There are {len(shortest_path)} shortest route detected")
            print(f"shortest route:", " -> ".join(shortest_path))
    
    # 绘图
    if plot == True:
        if full_route:
            count = 1
            print("...Saving Figures...")
            for m in tqdm(shortest_path):
                pos = nx.spring_layout(G)
                plt.figure(figsize=(20, 20))
                nx.draw_networkx_nodes(G, pos, node_size=30, node_color="#82B0D2", label=True, alpha=1)
                nx.draw_networkx_edges(G, pos, width=0.2, edge_color="#82B0D2", alpha=1)
                path_edges = [(m[i], m[i + 1]) for i in range(len(m) - 1)]
                path_nodes = m
                node_colors = ['#ec4347' if node in [start, end] else 'orange' for node in path_nodes]
                node_size = [400 if node in [start, end] else 300 for node in path_nodes]
                nx.draw_networkx_nodes(G, pos, nodelist=path_nodes, node_color=node_colors, node_size=node_size)
                shortest_path_labels = {node: node for node in path_nodes}
                nx.draw_networkx_labels(G, pos, labels=shortest_path_labels, font_size=7)
                nx.draw_networkx_edges(G, pos, edgelist=path_edges, width=1.8, edge_color='orange', arrows=True, arrowstyle='->')
                plt.axis('off')
                plt.savefig(f"{output}/path from {start} to {end} route {count + 1}.pdf")
                plt.close()
                count += 1
        else:
            pos = nx.spring_layout(G)
            plt.figure(figsize=(20, 20))
            nx.draw_networkx_nodes(G, pos, node_size=30, node_color="#82B0D2", label=True, alpha=1)
            nx.draw_networkx_edges(G, pos, width=0.2, edge_color="#82B0D2", alpha=1)
            path_edges = [(shortest_path[i], shortest_path[i + 1]) for i in range(len(shortest_path) - 1)]
            path_nodes = shortest_path
            node_colors = ['#ec4347' if node in [start, end] else 'orange' for node in path_nodes]
            node_size = [400 if node in [start, end] else 300 for node in path_nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=path_nodes, node_color=node_colors, node_size=node_size)
            shortest_path_labels = {node: node for node in path_nodes}
            nx.draw_networkx_labels(G, pos, labels=shortest_path_labels, font_size=7)
            nx.draw_networkx_edges(G, pos, edgelist=path_edges, width=1.8, edge_color='orange', arrows=True, arrowstyle='->')
            plt.axis('off')
            plt.savefig(f"{output}/restrict path from {start} to {end}.pdf")
            plt.close()
    else:
        pass


if __name__ == "__main__":
    '''
    record = open(f"data/short_path.txt", "a")
    for i in tqdm(range(1,404)):
        for j in range(1, 404):
            plot_graph_short_path('/Users/wangjingran/Desktop/PTEN mutation/pten_wt_edge.txt', 'Figure/path_short', str(i), str(j),record=record, plot=False)
    record.close()
    '''
    '''
    for i in range(1,404):
        plot_graph_short_path('/Users/wangjingran/Desktop/PTEN mutation/pten_wt_edge.txt', '/Users/wangjingran/Desktop/P_figure/path_short/101', '101', str(i), plot=True)
    '''
    graph_short_path('/Users/wangjingran/Desktop/WT_occupacy_network.ncol', '/Users/wangjingran/Desktop/P_figure/path_short', '88', '915',length = False, full_route = True)

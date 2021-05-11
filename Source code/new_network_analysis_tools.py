import networkx as nx
import random
import numpy as np
import matplotlib.pylab as plt
def read_network(network_file):
    """
    Reads a network from an external file.
    从外部文件读取网络

    * The edgelist must be provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2
    *edgelist必须以制表符分隔的形式提供。这个表的前两列将被解释为交互作用gene1<==>gene2

    * Lines that start with '#' will be ignored
    *以#号开头的行将被忽略
    """
    """创建一个空的无向图"""
    G = nx.Graph()
    for line in open(network_file, 'r'):
        # lines starting with '#' will be ignored
        if line[0] == '#':
            continue
        # The first two columns in the line will be interpreted as an
        # interaction gene1 <=> gene2
        # 这个表的前两列将被解释为交互作用gene1<==>gene2
        line_data = line.strip().split('\t')
        # line.strip()删除数据中的换行符 .split('\t')遇到四个空格就隔开
        node1 = line_data[0]
        node2 = line_data[1]
        # 定义node1为一组数据的前一项数据
        # 定义node2为一组数据的后一项数据
        G.add_edge(node1, node2)
        # 把数据作为一组加入到图中
    return G
    # 将数值返回给G 完成循环
def read_gene_list(gene_file):
    """
    Reads a list genes from an external file.

    * The genes must be provided as a table. If the table has more
    than one column, they must be tab-separated. The first column will
    be used only.

    * Lines that start with '#' will be ignored
    """

    genes_set = set()
    # 创建一个无序的不重复元素集，可进行关系测试，删除重复数据还可计算交集差集并集等
    for line in open(gene_file, 'r'):
        # lines starting with '#' will be ignored
        if line[0] == '#':
            continue
        # the first column in the line will be interpreted as a seed
        # gene:
        line_data = line.strip().split('\t')
        # line.strip()删除数据中的换行符.split('\t')遇到四个空格就隔开
        gene = line_data[0]
        # 只读取仅有的一号位置数据作为基因
        genes_set.add(gene)
        # 加入到set里面

    return genes_set
    # 返回genes_结束循环
# *********************************************************************
def get_closeness_centrality(network_file,gene_file):
    G = read_network(network_file)
    gene_set = read_gene_list(gene_file)
    g = nx.subgraph(G, gene_set)
    closeness_centrality_sum = 0
    closeness_centrality = list(nx.closeness_centrality(g).values())
    for i in range(0, len(closeness_centrality)):
        closeness_centrality_sum = closeness_centrality_sum + closeness_centrality[i]
    closeness_centrality_result = closeness_centrality_sum / len(closeness_centrality)
    return closeness_centrality_result
def get_random_closeness_centrality(network_file,random_times,gene_file):
    G = read_network(network_file)
    result_list=[]
    random_times = int(float(random_times))
    all_genes=G.nodes()
    gene_set = read_gene_list(gene_file)
    number_of_seed_genes = len(gene_set & set(all_genes))
    for i in range(0,random_times):
        random_seeds=set(np.random.choice(all_genes,number_of_seed_genes,replace=False))
        g1 = nx.subgraph(G,random_seeds)
        random_closeness_centrality_sum = 0
        random_closeness_centrality = list(nx.closeness_centrality(g1).values())
        for i in range(0, len(random_closeness_centrality)):
            random_closeness_centrality_sum = random_closeness_centrality_sum + random_closeness_centrality[i]
        random_closeness_centrality_result = random_closeness_centrality_sum / len(random_closeness_centrality)
        result_list.append(random_closeness_centrality_result)
    rand_closeness_centrality_result = 0
    for i in range(0,len(result_list)):
        rand_closeness_centrality_result = rand_closeness_centrality_result+result_list[i]
    rand_closeness_centrality_result = rand_closeness_centrality_result/len(result_list)
    closeness_centrality_result=get_closeness_centrality(network_file,gene_file)
    pvalue_list=[]
    for i in range(0,len(result_list)):
        if result_list[i]>closeness_centrality_result:
            pvalue_list.append(result_list[i])
    Pvalue = len(pvalue_list)/len(result_list)
    def all_list(arr):
        result2 = {}
        for i in set(arr):
            result2[i] = arr.count(i) / len(arr)
        return result2
    result = all_list(result_list)
    x = result.keys()
    y = result.values()
    plt.xlabel('The closeness centrality_percentage')
    plt.ylabel('closeness_centrality')
    plt.bar(x, y, width= 0.001, linewidth=0.005)
    plt.title('The histogram of the closeness_centrality_percentage')
    plt.savefig('The_closeness_centrality_final.pdf')
    plt.close()
    return rand_closeness_centrality_result,Pvalue
def get_Clustering_coefficient(network_file,gene_file):
    G = read_network(network_file)
    gene_set = read_gene_list(gene_file)
    g = nx.subgraph(G, gene_set)
    # print(nx.clustering(g))
    clustering_coefficient_sum = 0
    clustering_coefficient = list(nx.clustering(g).values())
    for i in range(0,len(clustering_coefficient)):
        clustering_coefficient_sum = clustering_coefficient_sum +clustering_coefficient[i]
    clustering_coefficient_result = clustering_coefficient_sum / len(clustering_coefficient)
    # print("clustering_coefficient的结果为：%s"%clustering_coefficient_result)
    return clustering_coefficient_result
def get_random_clustering_coefficient(network_file,gene_file,random_times):
    G = read_network(network_file)
    result_list = []
    random_times = int(float(random_times))
    all_genes = G.nodes()
    gene_set = read_gene_list(gene_file)
    number_of_seed_genes = len(gene_set & set(all_genes))
    for i in range(0, random_times):
        random_seeds = set(np.random.choice(all_genes, number_of_seed_genes, replace=False))
        g1 = nx.subgraph(G, random_seeds)
        random_clustering_coefficient_sum = 0
        random_clustering_coefficient = list(nx.clustering(g1).values())
        for i in range(0, len(random_clustering_coefficient)):
            random_clustering_coefficient_sum = random_clustering_coefficient_sum + random_clustering_coefficient[i]
        random_clustering_coefficient_result = random_clustering_coefficient_sum / len(random_clustering_coefficient)
        result_list.append(random_clustering_coefficient_result)
    rand_random_clustering_coefficient_result = 0
    for i in range(0, len(result_list)):
        rand_random_clustering_coefficient_result = rand_random_clustering_coefficient_result + result_list[i]
    rand_random_clustering_coefficient_result = rand_random_clustering_coefficient_result / len(result_list)
    clustering_coefficient_result = get_Clustering_coefficient(network_file, gene_file)
    pvalue_list = []
    for i in range(0, len(result_list)):
        if result_list[i] > clustering_coefficient_result:
            pvalue_list.append(result_list[i])
    Pvalue = len(pvalue_list) / len(result_list)
    def all_list(arr):
        result2 = {}
        for i in set(arr):
            result2[i] = arr.count(i) / len(arr)
        return result2

    result = all_list(result_list)
    x = result.keys()
    y = result.values()
    plt.xlabel('The clustering coefficient_percentage')
    plt.ylabel('clustering coefficient')
    plt.bar(x, y, width=0.001, linewidth=0.005)
    plt.title('The histogram of the clustering coefficient_percentage')
    plt.savefig('The_clustering_coefficient_final.pdf')
    plt.close()
    return rand_random_clustering_coefficient_result,Pvalue
def get_betweenness_centrality(network_file,gene_file):
    G = read_network(network_file)
    gene_set = read_gene_list(gene_file)
    g = nx.subgraph(G, gene_set)
    betweenness_centarlity_sum = 0
    betweenness_centarlity = list(nx.betweenness_centrality_source(g).values())
    for i in range(0, len(betweenness_centarlity)):
        betweenness_centarlity_sum = betweenness_centarlity_sum + betweenness_centarlity[i]
    betweenness_centarlity_result = betweenness_centarlity_sum / len(betweenness_centarlity)
    # print("clustering_coefficient的结果为：%s" % betweenness_centarlity_result)
    return betweenness_centarlity_result
def get_random_betweenness_centrality(network_file,gene_file,random_times):
    G = read_network(network_file)
    betw_result_list = []
    random_times = int(float(random_times))
    all_genes = G.nodes()
    gene_set = read_gene_list(gene_file)
    number_of_seed_genes = len(gene_set & set(all_genes))
    for i in range(0, random_times):
        random_seeds = set(np.random.choice(all_genes, number_of_seed_genes, replace=False))
        g1 = nx.subgraph(G, random_seeds)
        random_betweenness_centrality_sum = 0
        random_betweenness_centrality = list(nx.clustering(g1).values())
        for i in range(0, len(random_betweenness_centrality)):
            random_betweenness_centrality_sum = random_betweenness_centrality_sum + random_betweenness_centrality[i]
        random_betweenness_centrality_result = random_betweenness_centrality_sum / len(random_betweenness_centrality)
        betw_result_list.append(random_betweenness_centrality_result)
    rand_betweenness_centrality_result = 0
    for i in range(0, len(betw_result_list)):
        rand_betweenness_centrality_result = rand_betweenness_centrality_result + betw_result_list[i]
    rand_betweenness_centrality_result = rand_betweenness_centrality_result / len(betw_result_list)
    betweenness_centrality_result = get_betweenness_centrality(network_file, gene_file)
    pvalue_list = []
    for i in range(0, len(betw_result_list)):
        if betw_result_list[i] > betweenness_centrality_result:
            pvalue_list.append(betw_result_list[i])
    Pvalue = len(pvalue_list) / len(betw_result_list)
    def all_list(arr):
        result2 = {}
        for i in set(arr):
            result2[i] = arr.count(i) / len(arr)
        return result2

    result = all_list(betw_result_list)
    x = result.keys()
    y = result.values()
    plt.xlabel('The betweenness centrality_percentage')
    plt.ylabel('betweenness centrality')
    plt.bar(x, y, width=0.001, linewidth=0.005)
    plt.title('The histogram of the betweenness centrality_percentage')
    plt.savefig('The_betweenness_centrality_final.pdf')
    plt.close()
    return rand_betweenness_centrality_result,Pvalue
# *********************************************************************


#! /usr/bin/env python

"""
# -----------------------------------------------------------------------
#
# seperation.py
#
# by Joerg Menche
# Last Modified: 2014-12-06
#
# This code determines the network-based distance and sepration for
# two given sets of nodes on given network as described in
# 此代码确定给定网络上两个给定节点集的基于网络的距离和分离，如中所述
#
# Uncovering Disease-Disease Relationships Through The Human
# Interactome
# 通过人类相互作用体发现疾病与疾病的关系
#
# by Joerg Menche, Amitabh Sharma, Maksim Kitsak, Susan Dina
#    Ghiassian, Marc Vidal, Joseph Loscalzo & Albert-Laszlo Barabasi
#
#
# -----------------------------------------------------------------------
#
#
# This program will calculate the network-based distance d_AB and
# separation s_AB between two gene sets A and B.
# 该程序将计算两个基因集A和B之间基于网络的距离d_AB和分离s_AB。
#
# * Required input:（需求输入）
#
#   two files containing the gene sets A and B. The file must be in
#   form of a table, one gene per line. If the table contains several
#   columns, they must be tab-separated, only the first column will be
#   used. See the two files MS.txt and PD.txt for valid examples (they
#   contain genes for multiple sclerosis and peroxisomal disorders,
#   respectively).
#  包含基因集A和B的两个文件。该文件必须以表的形式存在，每行一个基因。
#  如果表包含多个列，它们必须用制表符分隔，只使用第一列。看这两个文件MS.txt
#  以及PD.txt文件作为有效示例，它们分别包含多发性硬化与过氧化物酶体疾病的有效基因。
# * Optional input:  （可选择输入）
#
#   - file containing an interaction network. If now file is given, the
#     default network \"interactome.tsv\" will be used instead. The file
#     must contain an edgelist provided as a tab-separated table. The
#     first two columns of the table will be interpreted as an
#     interaction gene1 <==> gene2
#  包含交互网络的文件。如果指定了now file，则默认网络\“Interactiome.tsv\“将被替代。
#  文件必须包含以制表符分隔的表形式提供的edgelist。表的前两列将被解释为交互gene1<=>gene2
#
#  - filename for the output. If none is given,
#    \"separation_results.txt\" will be used
#   对于输出的文件名，如果没有输入该参数，将会使用separetion_results.txt
#
# Here's an example that should work, provided the files are in the same
# directory as this python script:
# 下面是一个可以工作的示例，前提是文件与此python脚本位于同一目录中：
#
# ./separation.py -n interactome.tsv --g1 MS.txt --g2 PD.txt -o output.txt
#
#
# -----------------------------------------------------------------------
"""

""" networkx 用于构建和操作复杂的图结构，提供分析图的算法。"""
""" numpy 提供了python对多维数组对象的支持 """
""" optparse 它功能强大，而且易于使用，可以方便地生成标准的、符合Unix/Posix 规范的命令行说明"""
""" sys模块 提供了一系列有关Python运行环境的变量和函数。"""
import networkx as nx
import numpy as np
import optparse
import sys
import random
import matplotlib.pylab as plt
from matplotlib import mlab
from matplotlib import rcParams
from decimal import Decimal

"""
# =============================================================================

           S T A R T   D E F I N I T I O N S 

# =============================================================================
"""


# 已完成print_usage的注释
# =============================================================================
def print_usage(option, opt, value, parser):
    """定义一个帮助手册 """

    usage_message = """

# ----------------------------------------------------------------------

This program will calculate the network-based distance d_AB and
separation s_AB between two gene sets A and B.

* Required input:

  two files containing the gene sets A and B. The file must be in form
  of a table, one gene per line. If the table contains several
  columns, they must be tab-separated, only the first column will be
  used. See the two files MS.txt and PD.txt for valid examples (they
  contain genes for multiple sclerosis and peroxisomal disorders,
  respectively).

* Optional input:  

  - file containing an interaction network. If now file is given, the
    default network \"interactome.tsv\" will be used instead. The file
    must contain an edgelist provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2

 - filename for the output. If none is given,
   \"separation_results.txt\" will be used


Here's an example that should work, provided the files are in the same
directory as this python script:

./separation.py -n interactome.tsv --g1 MS.txt --g2 PD.txt -o output.txt

# ----------------------------------------------------------------------

    """

    print(usage_message)

    sys.exit()


# =============================================================================
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
    # 网络包括的节点数与相关联数
    return G
    # 将数值返回给G 完成循环


# =============================================================================
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
    # 打印输出基因数和文件名

    return genes_set
    # 返回genes_结束循环


# =============================================================================
def remove_self_links(G):
    # 该函数用于删除自我循环的数据例如 a-a
    sl = G.selfloop_edges()
    # 将会自我循环的数据定义给sl
    G.remove_edges_from(sl)
    # remove_edges_from删除sl


# =============================================================================
def get_pathlengths_for_single_set(G, given_gene_set):
    # 该对象作用为计算单个基因集在网络中的路径长度
    """
    calculate the shortest paths of a given set of genes in a
    given network. The results are stored in a dictionary of
    dictionaries:
    all_path_lenghts[gene1][gene2] = l
    with gene1 < gene2, so each pair is stored only once!
    计算给定网络中给定基因集的最短路径，结果储存在字典中
    PARAMETERS:（范围）
    -----------
        - G: network 基因网络
        - gene_set: gene set for which paths should be computed
        为其计算路径的基因集

    RETURNS:
    --------
        - all_path_lenghts[gene1][gene2] = l for all pairs of genes
          with gene1 < gene2（对于gene1<gene2的所有基因对）

    """

    # remove all nodes that are not in the network
    # 删除所有不在网络中的节点
    all_genes_in_network = set(G.nodes())
    # 将所有的nodes使用set()定义给all_genes_in_network
    gene_set = given_gene_set & all_genes_in_network
    # 使用&运算符将目的基因与所有的基因相同的定义给gene_set

    all_path_lenghts = {}

    # calculate the distance of all possible pairs
    # 计算所有可能对的距离
    for gene1 in gene_set:
        if not all_path_lenghts.__contains__(gene1):
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set:
            if gene1 < gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    all_path_lenghts[gene1][gene2] = l
                    # print(all_path_lenghts)
                except:
                    continue
    # 首先从目的基因（如MS.txt）中获取第一个基因，如果在基因网络中也有，则放在字典前如{gene1:{}}
    # 之后从总基因网络中获取与它相连的基因并计算而这最短路径 nx.shortest_path_length
    # 将结果方在字典后如{gene1:{gene2:length}}
    return all_path_lenghts


# =============================================================================
def get_pathlengths_for_single_set_of_random(G, given_gene_set, sims):
    """
        calculate the shortest paths of given set of genes from 1000 random samples
        in a given network. The results are stored in a dictionary of dictionaries:
        all_path_lenghts[gene1][gene2] = l
        with gene1 < gene2, so each pair is stored only once!
        计算给定网络中给定基因集的最短路径，结果储存在字典中
        PARAMETERS:（范围）
        -----------
            - G: network 基因网络
            - gene_set: gene set for which paths should be computed
            为其计算路径的基因集

        RETURNS:
        --------
            - all_path_lenghts[gene1][gene2] = l for all pairs of genes
              with gene1 < gene2（对于gene1<gene2的所有基因对）

        """
    all_genes_in_network = set(G.nodes())
    # 将所有的nodes使用set()定义给all_genes_in_network
    number = len(given_gene_set)
    # print(number)
    random_pathlengths_list = []
    for i in range(1, sims + 1):
        if i % 100 == 0:
            sys.stdout.write("> random simulation [%s of %s]\r" % (i, sims))
            sys.stdout.flush()
        rand_seeds = set(random.sample(all_genes_in_network, number))
        all_path_lenghts = get_pathlengths_for_single_set(G, rand_seeds)
        all_distances = []
        # going through all gene pairs 检查所有基因对
        for geneA in rand_seeds:
            all_distances_A = []
            for geneB in rand_seeds:

                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                # 我必须检查哪个基因“更小”，以便知道在全路径长度dict中从哪里查找这对基因的距离
                if geneA < geneB:
                    if all_path_lenghts[geneA].__contains__(geneB):
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                else:
                    if all_path_lenghts[geneB].__contains__(geneA):
                        all_distances_A.append(all_path_lenghts[geneB][geneA])

            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)

        # calculate mean shortest distance 计算平均最短距离 numpy包np
        mean_shortest_distance_fir = np.mean(all_distances)
        mean_shortest_distance = Decimal(mean_shortest_distance_fir).quantize(Decimal("0.00"))
        random_pathlengths_list.append(mean_shortest_distance_fir)
    return random_pathlengths_list


# =============================================================================
def get_pathlengths_for_two_sets(G, given_gene_set1, given_gene_set2):
    # 该对象作用是获取两组基因集在基因网络中的路径长度
    """
    calculate the shortest paths between two given set of genes in a
    given network. The results are stored in a dictionary of
    dictionaries: all_path_lenghts[gene1][gene2] = l with gene1 <
    gene2, so each pair is stored only once!
    计算给定网络中两组给定基因之间的最短路径。
    结果存储在字典字典中：
        all_path_lenghts[gene1][gene2]=l with  gene1<gene2
        每个对只存储一次！

    PARAMETERS:
    -----------
        - G: network
        - gene_set1/2: gene sets for which paths should be computed

    RETURNS:
    --------
        - all_path_lenghts[gene1][gene2] = l for all pairs of genes
          with gene1 < gene2

    """

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network

    all_path_lenghts = {}

    # calculate the distance of all possible pairs
    for gene1 in gene_set1:
        if not all_path_lenghts.__contains__(gene1):
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set2:
            if gene1 != gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    if gene1 < gene2:
                        all_path_lenghts[gene1][gene2] = l
                    else:
                        if not all_path_lenghts.__contains__(gene2):
                            all_path_lenghts[gene2] = {}
                        all_path_lenghts[gene2][gene1] = l
                except:
                    continue

    return all_path_lenghts


# =============================================================================
def calc_single_set_distance(G, given_gene_set):
    """
    Calculates the mean shortest distance for a set of genes on a
    given network    计算给定网络上一组基因的平均最短距离


    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which distance will be computed

    RETURNS:
    --------
         - mean shortest distance 平均最短距离

    """

    # remove all nodes that are not in the network, just to be safe
    # 为了安全起见，删除所有不在网络中的节点
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    # get the network distances for all gene pairs:获取所有基因对的网络距离：
    all_path_lenghts = get_pathlengths_for_single_set(G, gene_set)

    all_distances = []

    # going through all gene pairs 检查所有基因对
    for geneA in gene_set:

        all_distances_A = []
        for geneB in gene_set:

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            # 我必须检查哪个基因“更小”，以便知道在全路径长度dict中从哪里查找这对基因的距离
            if geneA < geneB:
                if all_path_lenghts[geneA].__contains__(geneB):
                    all_distances_A.append(all_path_lenghts[geneA][geneB])
            else:
                if all_path_lenghts[geneB].__contains__(geneA):
                    all_distances_A.append(all_path_lenghts[geneB][geneA])

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # calculate mean shortest distance 计算平均最短距离 numpy包np
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance

def calc_single_set_distance(G, given_gene_set):
    """
    Calculates the mean shortest distance for a set of genes on a
    given network    计算给定网络上一组基因的平均最短距离


    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which distance will be computed

    RETURNS:
    --------
         - mean shortest distance 平均最短距离

    """

    # remove all nodes that are not in the network, just to be safe
    # 为了安全起见，删除所有不在网络中的节点
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    # get the network distances for all gene pairs:获取所有基因对的网络距离：
    all_path_lenghts = get_pathlengths_for_single_set(G, gene_set)

    all_distances = []

    # going through all gene pairs 检查所有基因对
    for geneA in gene_set:

        all_distances_A = []
        for geneB in gene_set:

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            # 我必须检查哪个基因“更小”，以便知道在全路径长度dict中从哪里查找这对基因的距离
            if geneA < geneB:
                if all_path_lenghts[geneA].__contains__(geneB):
                    all_distances_A.append(all_path_lenghts[geneA][geneB])
            else:
                if all_path_lenghts[geneB].__contains__(geneA):
                    all_distances_A.append(all_path_lenghts[geneB][geneA])

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # calculate mean shortest distance 计算平均最短距离 numpy包np
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance
# =============================================================================
def calc_set_pair_distances(G, given_gene_set1, given_gene_set2):
    """
    Calculates the mean shortest distance between two sets of genes on
    a given network 计算给定网络上两组基因之间的平均最短距离

    PARAMETERS:
    -----------
        - G: network
        - gene_set1/2: gene sets for which distance will be computed

    RETURNS:
    --------
         - mean shortest distance

    """

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network

    # get the network distances for all gene pairs: 获取所有基因对的网络距离：
    all_path_lenghts = get_pathlengths_for_two_sets(G, gene_set1, gene_set2)

    all_distances = []

    # going through all pairs starting from set 1  从集合1开始检查所有对
    for geneA in gene_set1:

        all_distances_A = []
        for geneB in gene_set2:

            # the genes are the same, so their distance is 0
            if geneA == geneB:
                all_distances_A.append(0)

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass

                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass

        # 如果有多个distances_A 则求取其最小值并添加到all_distances
        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # going through all pairs starting from disease B 从疾病B开始检查所有配对
    for geneA in gene_set2:

        all_distances_A = []
        for geneB in gene_set1:

            # the genes are the same, so their distance is 0
            if geneA == geneB:
                all_distances_A.append(0)

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass

                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # calculate mean shortest distance 计算平均最短距离
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance


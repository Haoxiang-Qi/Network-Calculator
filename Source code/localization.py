#! /usr/bin/env python

"""
# -----------------------------------------------------------------------
#
# localization.py
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
# 通过人类相互作用组揭示疾病与疾病的关系
#
# by Joerg Menche, Amitabh Sharma, Maksim Kitsak, Susan Dina
#    Ghiassian, Marc Vidal, Joseph Loscalzo & Albert-Laszlo Barabasi
#
#
# -----------------------------------------------------------------------
#
#
# This program will calculate the size of the largest connected
# component S and mean shortest distance <d_s> for a given gene
# set. It will also compute the expected lcc size for the same number
# of randomly distributed genes.
# 此程序将计算最大连接的大小 给定基因组的组分S和平均最短距离<d_S>。它还将计算相同数量随机分布基因的预期lcc大小
# * Required input:
#
#   a file containing a gene set. The file must be in form of a table,
#   one gene per line. If the table contains several columns, they
#   must be tab-separated, only the first column will be used. See the
#   two files MS.txt and PD.txt for valid examples (they contain genes
#   for multiple sclerosis and peroxisomal disorders, respectively).
#
# * Optional input:
#
#   - file containing an interaction network. If now file is given, the
#     default network \"interactome.tsv\" will be used instead. The file
#     must contain an edgelist provided as a tab-separated table. The
#     first two columns of the table will be interpreted as an
#     interaction gene1 <==> gene2
#
#  - filename for the output. If none is given,
#    \"localiztion_results.txt\" will be used
#
#  - the number or random simulations can be chosen. Default is 1000,
#    which should run fast even for large gene sets and typically
#    gives good result.
#
# Here's an example that should work, provided the files are in the same
# directory as this python script:
#
# ./localization.py -n interactome.tsv -g PD.txt -o output.txt
#
#
# -----------------------------------------------------------------------
"""

import networkx as nx
import random
import numpy as np
import optparse
import sys
import matplotlib.pylab as plt
from collections import Counter
import separation as tools
from matplotlib import mlab
from matplotlib import rcParams
from random import choice
"""
# =============================================================================

           S T A R T   D E F I N I T I O N S 

# =============================================================================
"""


# 已完成print_usage注释.
# =============================================================================
def print_usage(option, opt, value, parser):
    usage_message = """

# ----------------------------------------------------------------------

This program will calculate the network-based localization for a given
gene set

* Required input:

  one files containing a gene set. The file must be in form of a
  table, one gene per line. If the table contains several columns,
  they must be tab-separated, only the first column will be used. See
  the two files MS.txt and PD.txt for valid examples (the contain
  genes for multiple sclerosis and peroxisomal disorders).

* Optional input:  

  - file containing an interaction network. If now file is given, the
    default network \"interactome.tsv\" will be used instead. The file
    must contain an edgelist provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2

  - filename for the output. If none is given,
    \"localiztion_results.txt\" will be used

  - the number or random simulations can be chosen. Default is 1000,
    which should run fast even for large gene sets and typically gives
    good result.


Here's an example that should work, provided the files are in the same
directory as this python script:

./localization.py -n interactome.tsv -g PD.txt -o output.txt

# ----------------------------------------------------------------------

    """

    print(usage_message)

    sys.exit()


# =================================================================================
def get_lcc_size(G, seed_nodes):
    """
    return the lcc size
    """
    # getting subgraph that only consists of the black_nodes 得到只包含节点的子图
    g = nx.subgraph(G, seed_nodes)
    if g.number_of_nodes() != 0:
        # get all components  获取所有组件
        components = list(sorted(nx.connected_components(g), key=len, reverse=True))
        return len(components[0])
    else:
        return 0


# =============================================================================
def get_edges_size(G, seed_nodes):
    """
    return the edges
    """
    # getting subgraph that only consists of the black_nodes 得到只包含节点的子图
    g = nx.subgraph(G, seed_nodes)
    return len(g.edges())


# =============================================================================
def get_random_mean_edges(G, gene_set, sims):
    # ============================================
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
    all_genes_in_network = set(G.nodes())
    tools.remove_self_links(G)
    # print(list(G))
    # =============================================
    all_genes = G.nodes()
    number_of_seed_genes = len(gene_set & set(all_genes))
    # print(number_of_seed_genes)
    # print(number_of_seed_genes)
    edges_list = []
    print("")
    for i in range(1, sims + 1):
        if i % 100 == 0:
            sys.stdout.write("> random simulation [%s of %s]\r" % (i, sims))
            sys.stdout.flush()
        rand_seeds = set(random.sample(all_genes,number_of_seed_genes))
        edges = get_edges_size(G, rand_seeds)
        edges_list.append(edges)

    # =============================================
    def all_list(arr):
        result2 = {}
        for i in set(arr):
            result2[i] = arr.count(i)

        return result2

    result2 = all_list(edges_list)
    # result2[result2.keys()]=result2.pop()
    # print(result)
    x = result2.keys()
    y = result2.values()
    plt.xlabel('edges')
    plt.ylabel('Frequency_percentage')
    plt.bar(x, y)
    plt.title('The histogram of the edges-Frequency_percentage')
    plt.savefig('edges_final.pdf')
    plt.close()
    # =============================================
    return edges_list


# =============================================================================
def get_random_comparison(G, gene_set, sims):
    # 获取随机比较

    """
    gets the random expectation for the lcc size for a given gene set
    by drawing the same number of genes at random from the network
    通过从网络中随机抽取相同数量的基因，获取给定基因集的lcc大小的随机期望值

    PARAMETERS:
    -----------
        - G       : network
        - gene_set: dito
        - sims    : number of random simulations

    RETURNS:
    --------
        - a string containing the results

    """
    # getting all genes in the network
    all_genes = G.nodes()
    number_of_seed_genes = len(gene_set & set(all_genes))  # 1000
    l_list = []

    # simulations with randomly distributed seed nodes 随机分布种子节点的模拟
    for i in range(1, sims + 1):
        # print out status  打印输出状态
        if i % 100 == 0:
            sys.stdout.write("> random simulation [%s of %s]\r" % (i, sims))
            sys.stdout.flush()
            # sys.stdout. 等价于print后加\n换行符

        # get random seeds
        rand_seeds = set(np.random.choice(all_genes, number_of_seed_genes,replace=False))

        # get rand lcc
        lcc = get_lcc_size(G, rand_seeds)
        l_list.append(lcc)

    # =================================================================
    # the histogram for the users
    def all_list(arr):
        result = {}
        for i in set(arr):
            result[i] = arr.count(i) / len(l_list)
        return result

    result = all_list(l_list)
    lcc = get_lcc_size(G, gene_set)
    x = result.keys()
    y = result.values()
    # print(y)
    plt.xlabel('lcc-Frequency_percentage')
    plt.ylabel('Frequency_percentage')
    plt.bar(x, y)
    plt.title('The histogram of the lcc-Frequency_percentage')
    plt.text(lcc - 4, 0.01, "lcc_size")

    plt.annotate("", xy=(lcc, 0.00002), xycoords='data',
                 xytext=(lcc, 0.01), textcoords='data',
                 arrowprops=dict(arrowstyle="->",
                                 connectionstyle="arc3"))

    plt.savefig('lcc_final.pdf')
    plt.close()

    # ==================================================================
    # get the actual value1
    lcc_observed = get_lcc_size(G, gene_set)

    # get the lcc z-score:
    l_mean = np.mean(l_list)
    l_std = np.std(l_list)

    if l_std == 0:
        z_score = 'not available'
    else:
        z_score = (1. * lcc_observed - l_mean) / l_std
    lcc = get_lcc_size(G, gene_set)
    new_list = []
    for i in l_list:
        if l_list[i] > lcc:
            new_list.append(l_list[i])
        i += 1
    P_value_of_lcc = (len(new_list)) / (len(l_list))
    new_edges_list = []
    # get_random_mean_edges(G, gene_set, sims)
    a = get_random_mean_edges(G, gene_set, sims)
    for i in a:
        if a[i] > edge_result:
            new_edges_list.append(a[i])
            i += 1
    P_value_of_edges = (len(new_edges_list) / len(a))
    results_message = """
> Random expectation:
> lcc [rand] = %s
> P-value of observed lcc = %s
> P-value of observed edges = %s
""" % (l_mean, P_value_of_lcc, P_value_of_edges)

    return results_message


def get_lcc_network(network_file, gene_set):
    dict1 = {}
    dict2 = {}
    for line in open(network_file, 'r'):
        # lines starting with '#' will be ignored
        if line[0] == '#':
            continue
        line_data = line.strip().split('\t')
        # line.strip()删除数据中的换行符 .split('\t')遇到四个空格就隔开
        node1 = line_data[0]
        node2 = line_data[1]
        dict1[node1] = node2
        dict2[node2] = node1
    g = nx.subgraph(G, gene_set)
    if g.number_of_nodes() != 0:
        # get all components  获取所有组件
        components = list(sorted(nx.connected_components(g), key=len, reverse=True))
        lcc_gene = list(components[0])
        lcc_network_forward = {}
        file = open('lcc_network.txt', 'a')
        print("The network of lcc:\n", file=file)
        for i in range(0, len(lcc_gene)):
            file = open('lcc_network.txt', 'a')
            print("%s--%s" % (lcc_gene[i], dict1[lcc_gene[i]]), file=file)
        for i in range(0, len(lcc_gene)):
            try:
                print("%s--%s" % (lcc_gene[i], dict2[lcc_gene[i]]), file=file)
                # print("matched [%d] "%i,file=file)
            except:
                # print("Matched over!",file=file)
                pass
        print("lcc networks 获取成功！", file=file)
        print("请将内容复制（CTRL+c）后从本文件夹中删除！", file=file)
        file.close()
def get_random_shortest_distance(G,gene_set):
    pass

"""
# =============================================================================

           E N D    O F    D E F I N I T I O N S 

# =============================================================================
"""

if __name__ == '__main__':

    # "Hey Ho, Let's go!" -- The Ramones (1976)

    # --------------------------------------------------------
    #
    # PARSING THE COMMAND LINE
    #
    # --------------------------------------------------------

    parser = optparse.OptionParser()

    parser.add_option('-u', '--usage',
                      help='print more info on how to use this script',
                      action="callback", callback=print_usage)

    parser.add_option('-n',
                      help='file containing the network edgelist [interactome.tsv]',
                      dest='network_file',
                      default='interactome.tsv',
                      type="string")

    parser.add_option('-g',
                      help='file containing gene set',
                      dest='gene_file',
                      default='none',
                      type="string")

    parser.add_option('-s',
                      help='number of random simulations [1000]',
                      dest='sims',
                      default='1000',
                      type="int")

    parser.add_option('-o',
                      help='file for results [separation_results.txt]',
                      dest='results_file',
                      default='localization_results.txt',
                      type="string")

    (opts, args) = parser.parse_args()

    network_file = opts.network_file
    gene_file = opts.gene_file
    results_file = opts.results_file
    sims = opts.sims

    # checking for input:
    if gene_file == 'none':
        error_message = """
        ERROR: you must specify an input file with a gene set, for example:
        ./localization.py -g MS.txt

        For more information, type
        ./localization.py --usage

        """
        print(error_message)
        sys.exit(0)

    if network_file == 'interactome.tsv':
        print('> default network from "interactome.tsv" will be used')

    # --------------------------------------------------------
    #
    # LOADING NETWORK and DISEASE GENES
    #
    # --------------------------------------------------------

    # read network
    G = tools.read_network(network_file)
    # get all genes ad remove self links
    all_genes_in_network = set(G.nodes())
    tools.remove_self_links(G)
    # read gene set
    gene_set_full = tools.read_gene_list(gene_file)
    # removing genes that are not in the network:
    gene_set = gene_set_full & all_genes_in_network
    if len(gene_set_full) != len(gene_set):
        print("> ignoring %s genes that are not in the network" % (
            len(gene_set_full - all_genes_in_network)))
        print("> remaining number of genes: %s" % (len(gene_set)))

    # --------------------------------------------------------
    #
    # CALCULATE NETWORK QUANTITIES
    #
    # --------------------------------------------------------

    # get lcc size S
    lcc = get_lcc_size(G, gene_set)

    print("\n> lcc size = %s" % (get_lcc_size(G, gene_set)))
    edge_result = get_edges_size(G, gene_set)
    mean = edge_result / len(gene_set)
    print("> edges = %s " % edge_result)
    print("> mean edges = %s" % mean)
    # get mean shortest distance
    d_s = tools.calc_single_set_distance(G, gene_set)
    print("> mean shortest distance = %s" % (d_s))

    results_message = """
> gene set from \"%s\": %s genes
> lcc size   S = %s
> diameter d_s = %s
""" % (gene_file, len(gene_set), lcc, d_s)

    # --------------------------------------------------------
    #
    # CALCULATE RANDOM COMPARISON
    #
    # --------------------------------------------------------

    results_message = results_message + get_random_comparison(G, gene_set, sims)
    print(results_message)
    get_lcc_network(network_file, gene_set)

    fp = open(results_file, 'w')
    fp.write(results_message)
    fp.close()
    print("> results have been saved to %s" % (results_file))












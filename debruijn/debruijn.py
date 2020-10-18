#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
import statistics
import argparse
import os
import sys
from random import randint
from random import randrange
import random
from operator import itemgetter
import networkx as nx
random.seed(9001)
#import matplotlib.pyplot as plt

__author__ = "Lynda"
__copyright__ = "Universite  de Paris "
__credits__ = ["Lynda"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lynda"
__email__ = "lyndamessad96@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
        Parameters:
        ---------
        path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file : str):
    """ This function is used to generate a sequence's generator from
the fastq file.
    Parameter:
    ---------
    fastq_file: str// fastq file with all sequences

    Return:
    ------
    sequences generator
    """
    with open(fastq_file) as file:
        for _ in enumerate(file): #parcourir les lines du fichier fastq
            yield next(file)[:-1]  # prend pas la derniere ligne du fichier
            next(file)
            next(file)


def cut_kmer(sequence : str, k_size : int):
    """
    This fiunction is used to genearte kmers from the sequences contained
    in the fastq file.

    Parameter:
    ---------
    sequence: str
    k_size: int // size of the kmer we want to split with

    Return:
    ------
    kmers's generator
    """
    for i in range(len(sequence)-k_size+1):
        yield sequence[i:i+k_size]


def build_kmer_dict(fastq_file : str, k_size : int):
    """
    This function is used to calculate  the number of occurence of
    kmers in a fq file.
    We use the functions cut_kmer() and read_fastq()

    Parameter:
    ---------
    fastq_file: str// fastq file with sequence
    k_size: int // size of the kmer we want to split with

    Return:
    ------
    kmer_dico: dictionary
        keys: kmer, values: nombre d'occurence of the kmer in the fastq file
    """
    dict_kmear = {}

    for seq in read_fastq(fastq_file):
        for k_mear in cut_kmer(seq, k_size):
            if k_mear not in dict_kmear.keys():
                dict_kmear[k_mear] = 1
            else:
                dict_kmear[k_mear] += 1

    return dict_kmear


def build_graph(k_mer_dict : dict):

    """
    This function creats the kmer's graph.
    It deppends of the kmer preffix,suffix and the weight.
    Weight: occurence of the kmer

    Parameter:
    ---------
    km_mear_dict: dictionary
        dictionary of kmer gets from build_kmer_dict() function

    Return:
    ------
    g: nx DiGraph : name = kmers_graph_kmer_size.png
    """
    graph = nx.DiGraph()

    for key in k_mer_dict:
        graph.add_edge(key[:-1], key[1:], weight = k_mer_dict[key])

    return graph


def get_starting_nodes(graph):
    """
    This function generats the starting nodes.

    Parameter:
    ---------
    graph: nx DiGraph// generated by the build_graph() function.

    Return:
    ------
    start_nodes: list // list of starting nodes found in the graph
    """
    start_nodes = []

    for node in graph.nodes:
        if len(graph.pred[node]) == 0:
        #Returns an iterator over predecessor nodes of n.
            start_nodes.append(node)

    return start_nodes

def get_sink_nodes(graph):
    """
    This function generats the sink nodes.

    Parameter:
    ---------
    graph: nx DiGraph// generated by the build_graph() function.

    Return:
    ------
    sink_nodes: list // list of sink nodes found in the graph
    """

    sink_nodes = []
    for node in graph.nodes:
        if len(graph.succ[node]) == 0:
            sink_nodes.append(node)

    return sink_nodes


def get_contigs(graph, start_nodes : list, sink_nodes : list):
    """
    This function generats a list of tuples that contains the
    contig and it's length.

    Parameter:
    ---------
    graph: nx DiGraph// generated by the build_graph() function.
    start_nodes : list// list of starting nodes found in the graph.
                        Generated by the function get_starting_nodes()

    sink_nodes: list// list of sink nodes found in the graph.
                    Generated by the function get_sink_nodes()

    Return:
    ------
    contigs_list: list of tuple// liste of tuple [(contig,len(contig))]
    """

    contigs_list = []
    for start in start_nodes: #start : source
        for sink in sink_nodes: #sink : target
            #Generate all simple paths in the graph G from
                # source(start) to target (skin)
            paths = list(nx.all_simple_paths(graph, start, sink))
            if paths:  #!= None ! -> returns True
                contig = paths[0][0]
                for i in range(1, len(paths[0])):
                    contig += paths[0][i][-1]
                contigs_list.append((contig, len(contig)))

    return contigs_list


def save_contigs(contigs : list, output_file : str):
    """
    This function creats the kmer's graph.
    It deppends of the kmer preffix,suffix and the weight.
    Weight: occurence of the kmer

    Parameter:
    ---------
    contig: list of tuple// liste of tuple [(contig,len(contig))]
    output_file : str // name of the output file

    Return:
    ------
    file : str // output file in fasta format generated by the functon fill()
    """
    with open(output_file, "w") as file:
        for contig, length in enumerate(contigs):
            line = ">contig_" + str(contig) + " len=" + str(length[1]) + "\n"
            file.write(line)
            file.write(fill(length[0]))
            file.write("\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std(liste : list):
    """
    This function calculats the standart deviation.
    If std > 0 so there is no variation between lenght and weight of paths.

    Parameter:
    ---------
    liste: list of values

    Return:
    ------
    float // the standart deviation of the liste
    """

    return statistics.stdev(liste)


def path_average_weight(graph, path : list):
    """
    This function calculats the average weight of the path gived in parameter.

    Parameter:
    ---------
    graph : object networkx DiGraph()
    path: list // path of de graph.

    Return:
    ------
    path with average weight
    """
    """weight = []

    for i in range(len(path) - 1):
        weight.append(graph.get_edge_data(path[i], path[i+1])["weight"])

    return statistics.mean(weight)

    """
    list_weight = []

    for u, v, e in graph.subgraph(path).edges(data=True):
        list_weight.append(e["weight"])

    return statistics.mean(list_weight)


def remove_paths(graph, paths : list, delete_entry_node : bool,
                 delete_sink_node : bool):
    """
    This function removes nodes (starting or sink) of path in graph.
    Parameter:
    ---------
    graph: object networkx DiGraph()// generated by the build_graph() function.
    path: list // path of de graph.
    delete_entry_node : boolean //if True, delete starting node of the path.
    delete_sink_node : boolean //if True, delete sink node of the path.

    Return:
    ------
    graph : object networkx DiGraph()// cleaned of unwanted paths.
    """
    for path in paths:
        for i in range(len(path)):

            if (i == 0 and not delete_entry_node) or \
                (i == len(path) - 1 and not delete_sink_node):
                continue
            graph.remove_node(path[i])

    return graph


def select_best_path(graph, paths : list, path_length : list, path_weight : list,
                     delete_entry_node=False, delete_sink_node=False):
    """
    This function select the best path for the graph.
    In this function the last 2 parametres are False (by default).

    Parameter:
    ---------
    graph : object networkx DiGraph()// generated by the build_graph() function.
    paths : list // path of de graph.
    path_length : list // indicates length of each path.
    path_weight : list // indicates average weight of each path.
    delete_entry_node : boolean // indicates if entry nodes will be deleted or not.
    delete_sink_node : boolean //  indicates if sink nodes will be deleted or not.

    Return:
    ------
    graph : object networkx DiGraph() // graph cleaned of unwanted paths.
    """
    best_path_len = 0
    best_path_index = -1 #out of range.

    for i in range(len(paths)):
        # Verify the biggest weight
        if  path_weight[i] == max(path_weight):
            best_path_len = path_length[i]
            best_path_index = i
            # compare by lenght if we have more than one best weight path
            if path_length[i] > best_path_len:
                best_path_len = path_length[i]
                best_path_index = i
            #in case the new best path have the same length than the previous
            # path than we do a random choice:
        elif path_length[i] == best_path_len:
            rand_path = randrange(0, 1, 1) #(start, stop, step)
            if rand_path: # True ->1
                best_path_len = path_length[i]
                best_path_index = i

    # Verify that we selected path, if not(just in case)we take a random choice
    if best_path_index == -1: #doesn't change (if True): no path selected.
        best_path_index = randint(0, len(paths))

    # Select the best path and delete the others
    removed_paths = paths[:best_path_index] + paths[best_path_index +1:]
    # remove from the graph
    graph = remove_paths(graph, removed_paths,
                          delete_entry_node,
                          delete_sink_node)

    return graph


def solve_bubble(graph, node_pred, node_succ):
    """
    This function is used to generate a graph cleaned of the bubble between
    the predecessor node and the successor node. Keep in bubble the best
    path selected.
    It uses the functions previously developed.

    Parameter:
    ---------
    graph: object networkx DiGraph()// generated by the build_graph() function.
    path_length : list // indicates length of each path.
    node_pred : predecessor node.
    node_succ : successor node.

    Return:
    ------
    graph : DiGraph // graph cleaned of bubble between node_succ and node_pred.
    """
    bubble_path_list = []
    path_length = [] #path length
    path_weight = [] #path weight

    for path in nx.all_simple_paths(graph, source=node_pred, target=node_succ):
        bubble_path_list.append(path)
        path_weight.append(path_average_weight(graph, path))
        path_length.append(len(path))

    # Keep the best path with select_best_path function
    return select_best_path(graph, bubble_path_list, path_length, path_weight)

def simplify_bubbles(graph):
    """
    This function simply the graph by taking off the bubbles.

    Parameter:
    ---------
    graph : object networkx DiGraph().

    Return:
    ------
    graph : object networkx DiGraph().
    """
    bubbles_nodes = []

    # Found out all bubbles
    for node in graph:
        node_ancesstors = [x for x in graph.predecessors(node)]
        if len(node_ancesstors) > 1: # Not empty
            ancestors = nx.lowest_common_ancestor(graph, node_ancesstors[0],
                        node_ancesstors[1])
            bubbles_nodes.append([ancestors, node])

    # Save the best path for each couple of bubbles.
    for nodes in bubbles_nodes:
        # Solve bubbles with previous function.
        graph = solve_bubble(graph, nodes[0], nodes[1])

    return graph


def solve_entry_tips(graph, entry_nodes):
    """
    This function is used to generate graph without unwanted entry nodes.
    This entry nodes must be gived as argument for the function.

    Parameter:
    ---------
    graph : object networkx DiGraph().
    entry_nodes : list of unwanted entry nodes.

    Return:
    ------
    graph : object networkx DiGraph().
    """
    paths_list, paths_length, paths_weight, descendants = [], [], [], []

    for node1 in entry_nodes:
        for node2 in nx.descendants(graph, node1):
            if len(graph.pred[node2]) >= 2: #at least 2 pred
                if node2 not in descendants:
                    descendants.append(node2)
                    #list of descendants nodes

    for node1 in entry_nodes:
        for node2 in descendants:
            #selcet simple path function
            for simple_path in nx.all_simple_paths(graph, node1, node2):
                #Get arguments to select best paths with the //
                # the select_best_path function
                paths_list.append(simple_path)
                paths_length.append(len(simple_path))
                paths_weight.append(path_average_weight(graph, simple_path))
        #Run function select_best_path avec delete_entry_node = True
        # because we solve entry tips.
        graph = select_best_path(graph, paths_list,
                                 paths_length, paths_weight,
                                delete_entry_node=True,
                                delete_sink_node=False)

    return graph

def solve_out_tips(graph, sink_nodes):
    """
    This function is used to generate graph without unwanted sink nodes.
    This output nodes must be gived as seconde argument for the function.

    Parameter:
    ---------
    graph : object networkx DiGraph().
    sink_nodes : list of unwanted entry nodes.

    Return:
    ------
    graph : object networkx DiGraph().
    """
    paths_list, paths_length, paths_weight, descendants = [], [], [], []

    for node1 in sink_nodes:
        for node2 in nx.ancestors(graph, node1):
            if len(graph.succ[node2]) >= 2 and node2 not in descendants:
                descendants.append(node2)

    for node in sink_nodes:
        for node2 in descendants:
            for simple_path in nx.all_simple_paths(graph, node2, node):
                #Get arguments to select best paths with the //
                # the select_best_path function

                paths_list.append(simple_path)
                paths_length.append(len(simple_path))
                paths_weight.append(path_average_weight(graph, simple_path))

        #Run function select_best_path avec delete_entry_node = True
        # because we solve entry tips.
        graph = select_best_path(graph, paths_list, paths_length, paths_weight,
                             delete_entry_node=False, delete_sink_node=True)

    return graph

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # 1. Bruijn's graph conception
        # 1.a kmer occuracy disctionary
    dict_kmer_occur = build_kmer_dict(args.fastq_file,args.kmer_size)

        # 1.b Buijn'tree conception
    graph = build_graph(dict_kmer_occur)

    # 2. Manipulation of Bruijn's graph
        #Contigs:
    contigs = get_contigs(graph,
                          get_starting_nodes(graph), get_sink_nodes(graph))

        # Save contigs in file
    save_contigs(contigs, args.output_file)

    graph = simplify_bubbles(graph)
if __name__ == '__main__':
    main()

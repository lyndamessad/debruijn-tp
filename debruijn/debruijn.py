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

import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Lynda"
__copyright__ = "Universite  de Paris "
__credits__ = ["Lynda"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "lyndamessad96@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
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

def read_fastq(fastq_file):
    """ This function is used to generate a sequence generator
        input: fastq file
        Return: sequence generator
    """
    with open(fastq_file) as file:
        for line in enumerate(file): #parcourir les lines du fichier fastq
            yield next(file)[:-1]  # On ne prend pas la derniere lignes du fichier 
            next(file)  # on saute les 2 lignes pour ne récupérer que les séquence 
            next(file) 

def cut_kmer(sequence, k_size):
    """ This fiunction is used to genearte kmears from the sequences.
    Input: sequence ( string) and k-mer (integer)
    return: k-mer generator
    """
    k = k_size
    for i in range(len(sequence)-k+1):
        yield (sequence[i:i+k])



def build_kmer_dic(fastq_file, k_size):
    """
    This function is used to calculate  the number of occurence of kmears in a fq file. 
    We use the functions cut_kmer() and read_fastq()
    """
    dict_kmear = {}
    k = k_size

    for seq in read_fastq(fastq_file): 
        for k_mear in cut_kmer(seq, k):
            if k_mear not in dict_kmear: 
           
                dict_kmear[k_mear] = 1
            else: 
                dict_kmear[k_mear] += 1

    return dict_kmear

def build_graph(k_mear_dict):
    """
    trouver les prefixes et suffixes pour le graph 
    """
    

    nodes = set() #pas de redondance pas besoin de boucle pr verifier 
    edges = {} #edge : weight 
    kmears = list(k_mear_dict)
    k = len(kmears[0])
    for kmear in kmears:
        suffix = kmear[1:]
        prefixe = kmear[:-1]
        nodes.add(prefixe)
        nodes.add(suffix)
        edges[(prefixe, suffix)] =  k_mear_dict[kmear]
    """
    print(edges.keys())
    G.add_nodes_from(nodes)
    G.add_edges_from(edges.keys(), weight = edges.values())

    nx.draw(G)
    plt.savefig("path.png")
    plt.show()
    """
    G=nx.DiGraph()
    G.add_edges_from(edges)

    
    
    print(G.nodes())
    nx.draw(G)
    plt.savefig("path.png")
    plt.show()
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    fastq_file=  vars(args)['fastq_file']
    sequences = read_fastq(fastq_file)
    dict_kmear_occur = build_kmer_dic(fastq_file,3)
    #print(dict_kmear_occur,"\n")
    graph = build_graph(dict_kmear_occur)

    #print(graph)

    #nodes, edges=  build_graph(dict_kmear_occur)
    #print(edges)


if __name__ == '__main__':
    main()


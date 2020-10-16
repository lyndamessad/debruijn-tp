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
__maintainer__ = "Lynda "
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
    """ This function is used to generate a sequence's generator from the fastq file. 
    Parameter: 
    ---------
    fastq_file: str// fastq file with all sequences 

    Return:
    ------
    sequences generator 
    """
    with open(fastq_file) as file:
        for line in enumerate(file): #parcourir les lines du fichier fastq
            yield next(file)[:-1]  # On ne prend pas la derniere lignes du fichier 
            next(file)  # on saute les 2 lignes pour ne récupérer que les séquence 
            next(file) 

def cut_kmer(sequence, k_size):
    """ 
    Input: sequence ( string) and k-mer (integer)
    return: k-mer generator
    """

    """
    This fiunction is used to genearte kmers from the sequences contained in the fastq file.

    Parameter: 
    ---------
    sequence: str 
    k_size: int // size of the kmer we want to split with 

    Return:
    ------
    kmers's generator
    """


    k = k_size
    for i in range(len(sequence)-k+1):
        yield (sequence[i:i+k])



def build_kmer_dic(fastq_file, k_size):
    """
    This function is used to calculate  the number of occurence of kmears in a fq file. 
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
    k = k_size

    for seq in read_fastq(fastq_file): 
        for k_mear in cut_kmer(seq, k):
            if k_mear not in dict_kmear.keys(): 
           
                dict_kmear[k_mear] = 0
            else: 
                dict_kmear[k_mear] += 1

    return dict_kmear


def build_graph(k_mear_dict):

    """
    This function creats the kmer's graph.
    It deppends of the kmer preffix,suffix and the weight. 
    Weight: occurence of the kmer 

    Parameter: 
    ---------
    km_mear_dict: dictionary
        dictionary of kmer gets from build_kmer_dic() function
    
    Return:
    ------
    graph: nx DiGraph : name = kmers_graph_kmer_size.png
    """

    k = vars(get_arguments())['kmer_size'] #kmer's size
    G=nx.DiGraph()
    for kmear, w in k_mear_dict.items():
        G.add_edge(kmear[:-1], kmear[1:], weight = w)
    plt.subplot(111)
    
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.savefig("kmers_graph_{}.png".format(k))  #graph name 
    
    return G
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
    k = vars(args)['kmer_size']
    sequences = read_fastq(fastq_file)
    dict_kmear_occur = build_kmer_dic(fastq_file,k)
    graph = build_graph(dict_kmear_occur)



if __name__ == '__main__':
    main()

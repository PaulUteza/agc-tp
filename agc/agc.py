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

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "UTEZA Paul"
__copyright__ = "CY Tech"
__credits__ = ["UTEZA Paul"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "UTEZA Paul"
__email__ = "utezapaul@eisti.eu"
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
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    Read amplicon file and generate FASTA sequences superior of a certain length
    :param amplicon_file: File with FASTA sequences
    :param minseqlen: minimum length of sequences we will consider
    """
    with gzip.open(amplicon_file, "rt") as file:
        sequence = ""
        sequence_list = []
        for line in file:
            if not line.lstrip().startswith('>'):
                sequence += line.strip()

            else:
                if len(sequence) != 0:
                    sequence_list.append(sequence)
                sequence = ""
    if len(sequence) != 0:
        sequence_list.append(sequence)
    for sequence in sequence_list:
        if len(sequence) >= minseqlen:
            yield sequence

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Generate unique sequences and corresponding occurence count with more occurences than mincount
    :param amplicon_file: File containing FASTA sequences
    :param minseqlen: minimum length of sequences we will consider
    :param mincount: minimum occurrence count of sequences we will consider
    """
    sequences = read_fasta(amplicon_file, minseqlen)
    sequence_list = list(sequences)
    seq_counter = Counter(sequence_list)
    seq_order = sorted(seq_counter, key=seq_counter.get, reverse=True)
    for seq in seq_order:
        if seq_counter[seq] >= mincount:
            yield [seq, seq_counter[seq]]


def get_chunks(sequence, chunk_size):
    """
    Create list of sub sequences of given size with a minimum of 4 sequences
    :param sequence: Sequence to split
    :param chunk_size: Size of sub sequences
    :return: Sub sequences of given size
    """
    if len(sequence)//chunk_size < 4:
        print('Not enough chunks')
        raise ValueError
    return [sequence[i:i + chunk_size] for i in range(0, len(sequence), chunk_size)
            if len(sequence[i:i + chunk_size]) >= chunk_size]


def get_unique(ids):
    """
    Return unique keys
    :param ids: List
    :return: unique keys
    """
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    """
    Return list of common elements from two lists
    :param lst1: First list to compare
    :param lst2: Second list to compare
    :return: List of unique common elements from both lists
    """
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
    Take a sequence and generate kmer cut in given size
    :param read: sequence
    :param kmer_size: size the kmer needs to be cut
    """
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:kmer_size + i]

def get_identity(alignment_list):
    """
    Return percentage of identity (similarity) from a list of 2 alignements
    :param alignment_list: List of alignements
    :return: Percentage of identity
    """
    count = 0
    for i in range(0, len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            count +=1
    return round(count/len(alignment_list[0])*100, 2)

def detect_chimera(perc_indentity_matrix):
    """
    Detect a chimera based on rules from an identity matrix
    :param perc_indentity_matrix: Identity matrix in percentage
    :return: Boolean depending if it's a chimera or not
    """
    seq_1 = []
    seq_2 = []
    std_perc = 0
    for perc_id in perc_indentity_matrix:
        std_perc += statistics.stdev(perc_id)
        seq_1.append(perc_id[0])
        seq_2.append(perc_id[1])

    if len(set(seq_1)) >= 2 or len(set(seq_2)) >= 2:
        if std_perc/len(perc_id) > 5.0:
            return True
    return False


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    Create dictionary of kmer with their sequence id as value
    :param kmer_dict: dictionary of kmer
    :param sequence: sequence to consider
    :param id_seq: id of the sequence
    :param kmer_size: size the kmer needs to be cut
    :return: Dictionary of kmer with their sequence id as value
    """
    cut_kmers = cut_kmer(sequence, kmer_size)
    for kmer in cut_kmers:
        if kmer in kmer_dict.keys():
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
    Find 8 most similar sequences of given sequence
    :param kmer_dict: dictionary of kmer
    :param sequence: sequence to consider
    :param kmer_size: size the kmer needs to be cut
    :return: 8 most similar sequences of given sequence
    """
    cut_kmers = cut_kmer(sequence, kmer_size)
    mates = []
    for kmer in cut_kmers:
        if kmer in kmer_dict.keys():
            mates.append(kmer_dict[kmer])
        else:
            continue
    mates = [item for sublist in mates for item in sublist]
    counter_mates = Counter(mates)
    best_mates_counter = counter_mates.most_common(8)
    best_mates = [x[0] for x in best_mates_counter]
    return best_mates


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Generate non chimera sequences
    :param amplicon_file: File containing FASTA sequences
    :param minseqlen: minimum length of sequences we will consider
    :param mincount: minimum occurence count of sequences we will consider
    :param chunk_size: Size of sub sequences
    :param kmer_size: size the kmer needs to be cut
    """
    dereplications = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    kmer_dict = {}
    no_chimera = []
    for seq_id, (dereplication, dereplication_count) in enumerate(dereplications):
        chunks = get_chunks(dereplication, chunk_size)
        mates = [search_mates(kmer_dict, chunk, kmer_size) for chunk in chunks]
        parents = mates[0]
        for mate in mates:
            parents = common(parents, mate)[:2]
        perc_identity_matrix = []
        if len(parents) == 2:
            perc_identity_matrix = [[] for _ in chunks]
            ref_chunks_1 = get_chunks(no_chimera[parents[0]], chunk_size)
            ref_chunks_2 = get_chunks(no_chimera[parents[1]], chunk_size)
            for idx, ref in enumerate(ref_chunks_1):
                perc_identity_matrix[idx] = get_identity(nw.global_align(chunks[idx], ref))
            for idx, ref in enumerate(ref_chunks_2):
                perc_identity_matrix[idx] = get_identity(nw.global_align(chunks[idx], ref))
        if not detect_chimera(perc_identity_matrix):
            no_chimera.append(dereplication)
            kmer_dict = get_unique_kmer(kmer_dict, dereplication, seq_id, kmer_size)
            yield [dereplication, dereplication_count]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Greedy clustering using non chimera sequence generator to find
    if another sequence is similar at 97%+ and more common
    :param amplicon_file: File containing FASTA sequences
    :param minseqlen: minimum length of sequences we will consider
    :param mincount: minimum occurence count of sequences we will consider
    :param chunk_size: Size of sub sequences
    :param kmer_size: size the kmer needs to be cut
    :return: List of OTU
    """
    pass

def fill(text, width=80):
    """
    Create correct FASTA format string
    :param text: Text to split in FASTA format
    :param width: max width per line
    :return: String in correct format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """
    Save output in file with FASTA format
    :param OTU_list: List of OTU sequences
    :param output_file: file to save
    """
    pass


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get args
    args = get_arguments()

    # Get file and variables
    file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size
    output_file = args.output_file



if __name__ == '__main__':
    main()

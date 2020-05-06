from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))

K = 25


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None

def create_kmer_from_read(read, k):
    kmers = []
    for i in range(len(read) - (k - 1)):
        kmers.append(read[i:(i + k)])
    return kmers

def create_kmers(reads, k):
    kmers = []
    for read in reads:
        kmers += create_kmer_from_read(read, k)
    return kmers

def expand_dict_to_list(dict, valid_keys = None):
    if valid_keys == None:
        valid_keys = dict.keys()
    result = []
    for key in valid_keys:
        for i in range(dict[key]):
            result.append(key)
    return result

def reduce_to_number_edges(kmers_to_count, valid_kmers = None):
    if valid_kmers == None:
        valid_kmers = kmers_to_count.keys()
    kmers_reduced = {}
    for kmer in valid_kmers:
        kmers_reduced[kmer] = round(kmers_to_count[kmer]/K)
    return kmers_reduced

def remove_errors(kmers):
    kmers_to_count = {}
    for kmer in kmers:
        if kmer in kmers_to_count:
            kmers_to_count[kmer] += 1
        else:
            kmers_to_count[kmer] = 1

    valid_kmers = [kmer for kmer in kmers_to_count.keys() if kmers_to_count[kmer] < K/5]
    valid_kmers = reduce_to_number_edges(kmers_to_count, valid_kmers)
    return expand_dict_to_list(kmers_to_count, valid_kmers)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    single_reads = [read for read_pair in input_reads for read in read_pair]

    kmers = remove_errors(create_kmers(single_reads[:5], K))

    # contigs = ['GCTGACTAGCTAGCTACGATCGATCGATCGATCGATCGATGACTAGCTAGCTAGCGCTGACT']

    contigs = kmers

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

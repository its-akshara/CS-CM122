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

def get_cycles(dic, cycles, start):
	if dic == {}:
		return {}, cycles, start
	curr = start
	c = []
	next = None

	while True:
		c.append(curr)
		if next not in out_deg.keys() or next not in in_deg.keys():	
			if curr not in dic.keys():
				return dic, [], start

			next = dic[curr].pop(0)
			curr = next
		else:
			break

	cycles.append(c)

	to_del = []
	for k in dic.keys():
		if dic[k]  == []:
			to_del.append(k)
	for k in to_del:
		del dic[k]

	if dic == {}:
		return dic, cycles, start

	for i in c:
		if i in dic.keys():
			dic, cycles, _ = get_cycles(dic, cycles, i)
			if dic == {}:
				return dic, cycles, start

	return dic, cycles, start

def get_maximal_path(path, cycle, rest):
	if cycle == []:
		return path, [], []

	for c in cycle:
		for r in rest:
			if c == r[0]:
				rest.remove(r)
				path, _, rest = get_maximal_path(path, r, rest)
		else:
			path.append(c)

	return path[:-1], [], rest

def create_debruijn(kmers):
    print(len(kmers))
    debruijn = {}
    for kmer in kmers:
        if kmer[:len(kmer)-1] not in debruijn:
            debruijn[kmer[:len(kmer)-1]] = [kmer[1:]]
        else:
            debruijn[kmer[:len(kmer)-1]].append(kmer[1:])

    max = 0
    for key in debruijn:
        if len(debruijn[key]) > max:
            max = len(debruijn[key])
    print("max={}".format(max))
    return debruijn

def generate_contigs(dic):    
    global in_deg
    global out_deg
    in_deg = {}
    out_deg = {}

    for k in dic.keys():
        out_deg[k] = len(dic[k])
        for v in dic[k]:
            if v not in in_deg.keys():
                in_deg[v] = 1
            else:
                in_deg[v] += 1
        
    for k in dic.keys():
        if k not in out_deg.keys():
            out_deg[k] = 0
        if k not in in_deg.keys():
            in_deg[k] = 0

        for v in dic[k]: 
            if v not in out_deg.keys():
                out_deg[v] = 0
            if v not in in_deg.keys():
                in_deg[v] = 0

    remove = []
    for k in in_deg.keys():
        if in_deg[k] == out_deg[k] == 1:
            remove.append(k)

    for k in remove:
        del in_deg[k]
        del out_deg[k]

    paths = []

    start = None
    for k in in_deg.keys():
        start = k   
        cycles = []
        dic, cycles, _ = get_cycles(dic, cycles, start)
        
        if cycles == []:
            continue

        for c in cycles:
            paths.append(c)
    
    return create_contigs_from_paths(paths)

def create_contigs_from_paths(paths):
    contigs = []
    for path in paths:
        ordered = path[0]
        for i in range(1,len(path)):
            ordered += path[i][-1]
        contigs.append(ordered)
    return contigs

def expand_dict_to_list(dict):
    result = []
    for key in dict:
        for i in range(dict[key]):
            result.append(key)
    return result

def reduce_to_number_edges(kmers_to_count):
    kmers_reduced = {}
    for kmer in kmers_to_count:
        # if kmers_to_count[kmer] >= 30 and kmers_to_count[kmer] < 60:
        #     kmers_reduced[kmer] = 1
        # elif kmers_to_count[kmer] < 90:
        #     kmers_reduced[kmer] = 2
        # else:
        #     print("A lot: {}".format(kmers_to_count[kmer]))
        #     kmers_reduced[kmer] = 3
        # occurences = int(kmers_to_count[kmer]/30)
        # if occurences > 0:
            # kmers_reduced[kmer] = occurences
        if kmers_to_count[kmer] > 3:
            kmers_reduced[kmer] = 1
    print("reduced len = {}".format(len(kmers_reduced)))
    return kmers_reduced

def create_kmers_to_count(kmers):
    kmers_to_count = {}
    for kmer in kmers:
        if kmer in kmers_to_count:
            kmers_to_count[kmer] += 1
        else:
            kmers_to_count[kmer] = 1
    print("count len = {}".format(len(kmers_to_count)))
    return kmers_to_count

def remove_errors(kmers):
    kmers_to_count = create_kmers_to_count(kmers)

    valid_kmers = reduce_to_number_edges(kmers_to_count)
    return expand_dict_to_list(valid_kmers)

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

    kmers = remove_errors(create_kmers(single_reads, K))
    debruijn = create_debruijn(kmers)
    contigs = generate_contigs(debruijn)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

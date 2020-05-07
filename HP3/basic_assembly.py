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

def delete_dead_debruijn_nodes(debruijn):
    delete_queue = []
    for node in debruijn:
        if len(debruijn[node]) == 0:
            delete_queue.append(node)
    for node in delete_queue:
        del debruijn[node]
    return debruijn

def get_cycles(debruijn_left, cycles, start):
	if debruijn_left == {}:
		return {}, cycles, start
	curr = start
	cycle = []
	next = None

	while True:
		cycle.append(curr)
		if next not in out_deg or next not in in_deg:	
			if curr not in debruijn_left:
				return debruijn_left, [], start

			next = debruijn_left[curr].pop(0)
			curr = next
		else:
			break
	cycles.append(cycle)

	debruijn_left = delete_dead_debruijn_nodes(debruijn_left)

	if debruijn_left == {}:
		return debruijn_left, cycles, start

	for node in cycle:
		if node in debruijn_left:
			debruijn_left, cycles, _ = get_cycles(debruijn_left, cycles, node)
			if debruijn_left == {}:
				return debruijn_left, cycles, start

	return debruijn_left, cycles, start

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
    debruijn = {}
    for kmer in kmers:
        if kmer[:len(kmer)-1] not in debruijn:
            debruijn[kmer[:len(kmer)-1]] = [kmer[1:]]
        else:
            debruijn[kmer[:len(kmer)-1]].append(kmer[1:])
    return debruijn

def calculate_degrees_debruijn(debruijn):
    global in_deg
    global out_deg
    in_deg = {}
    out_deg = {}

    for node in debruijn:
        out_deg[node] = len(debruijn[node])
        for neighbor in debruijn[node]:
            if neighbor not in in_deg:
                in_deg[neighbor] = 1
            else:
                in_deg[neighbor] += 1
        
    for node in debruijn:
        if node not in out_deg:
            out_deg[node] = 0
        if node not in in_deg:
            in_deg[node] = 0

        for neighbor in debruijn[node]: 
            if neighbor not in out_deg:
                out_deg[neighbor] = 0
            if neighbor not in in_deg:
                in_deg[neighbor] = 0

def generate_contigs(debruijn):    
    calculate_degrees_debruijn(debruijn)

    remove = []
    for k in in_deg:
        if in_deg[k] == out_deg[k] == 1:
            remove.append(k)

    for k in remove:
        del in_deg[k]
        del out_deg[k]

    paths = []

    for start_node in in_deg: 
        cycles = []
        debruijn, cycles, _ = get_cycles(debruijn, cycles, start_node)

        for cycle in cycles:
            paths.append(cycle)
    
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
        if kmers_to_count[kmer] > 3:
            kmers_reduced[kmer] = 1
    return kmers_reduced

def create_kmers_to_count(kmers):
    kmers_to_count = {}
    for kmer in kmers:
        if kmer in kmers_to_count:
            kmers_to_count[kmer] += 1
        else:
            kmers_to_count[kmer] = 1
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

import sys
import argparse
import numpy as np
import time
import zipfile

L = 30
MISMATCHES_ALLOWED = 2 # Number of mismatches allowed
CONSENSUS_MAJORITY = 175 # Number to get consensus that SNP located there

def create_subsequence_lookup(genome):
    subseq_to_index = {}

    for i in range(int(len(genome) - L/3 + 1)):
        seq = genome[i:int(i+L/3)]
        if seq in subseq_to_index:
            subseq_to_index[seq].append(i)
        else:
            subseq_to_index[seq] = [i]

    return subseq_to_index

def split_into_3(read):
    return ([read[i:int(i + L/3)] for i in range(0, len(read), int(L/3))])

def ref_start_pos(index, which_third):
    if which_third == 0:
        return index
    elif which_third == 1:
        return int(index - L/3)
    else:
        return int(index - (2*L)/3)

def find_pos_differences(ref, read, start_pos):
    diff = []

    for i in range(len(read)):
        if read[i] != ref[i]:
            diff.append([ref[i], read[i] , start_pos+i])

    return diff

def evaluates_indices(indices, read, ref, which_third):
    snps = []

    for i in range(len(indices)):
        ref_start = ref_start_pos(indices[i], which_third)
        ref_subseq = ref[ref_start:(ref_start+L)]

        if (len(ref_subseq)) == L:
            diff = find_pos_differences(ref_subseq, read, ref_start) 

            if len(diff) <= MISMATCHES_ALLOWED and len(diff) > 0:
                snps += diff

    return snps

# returns list of [OriginalAllele,SNP,Position]
def find_possible_snp_in_read(read, lookup, ref):
    possible_snps = []
    thirds = split_into_3(read)
    which_third = 0
    
    for third in thirds:
        if third in lookup:
            possible_indices = lookup[third]
            snps_for_third = (evaluates_indices(possible_indices, read, ref, which_third))
            possible_snps += snps_for_third

        which_third += 1
    
    return possible_snps

def count_occurences_possible_snps(snps):
    snp_pos_to_count = {}

    for snp_tuple in snps:
        if tuple(snp_tuple) in snp_pos_to_count:
            snp_pos_to_count[tuple(snp_tuple)] += 1
        else:
            snp_pos_to_count[tuple(snp_tuple)] = 1

    return snp_pos_to_count

def choose_majority(possibilities_to_count, majority):
    diffs = []
    # avg = 0
    for diff_tuple in possibilities_to_count:
        if possibilities_to_count[diff_tuple] >= majority:
            # print(possibilities_to_count[diff_tuple])
            diffs.append(list(diff_tuple))
            # avg += possibilities_to_count[diff_tuple]
    # print("AVG={}".format(avg/len(diffs)))
    return diffs

def remove_pos_duplicates(sorted_snps, snp_possibilities_to_count):
    pos_to_snp = {}
    final_snps = []
    for snp in sorted_snps:
        if snp[2] in pos_to_snp:
            prev_snp = (pos_to_snp[snp[2]][0], pos_to_snp[snp[2]][1], snp[2])
            prev_count = snp_possibilities_to_count[prev_snp]
            if prev_count < snp_possibilities_to_count[tuple(snp)]:
                final_snps.pop()
            else:
                continue
        pos_to_snp[snp[2]] = (snp[0], snp[1])
        final_snps.append(snp)
    return final_snps

def find_snps(reads, lookup, ref):
    possible_snps = []
    snps = []

    for read in reads:
        possible_snps += find_possible_snp_in_read(read, lookup, ref)

    snp_possibilities_to_count = count_occurences_possible_snps(possible_snps)
    snps = choose_majority(snp_possibilities_to_count, CONSENSUS_MAJORITY)
    snps.sort(key=lambda x:x[2])

    return remove_pos_duplicates(snps, snp_possibilities_to_count)

def kmer_comp(read):
    kmers = []
    for i in range(int(len(read) - L + 1)):
        kmers.append(read[i:int(i+L)])
    return kmers

def enumerate_reads(reads):
    kmers = []
    for read in reads:
        kmers += kmer_comp(read)
    return kmers

def convert_pairs_to_reads(paired_reads):
    return [read for read_pair in paired_reads for read in read_pair]

def is_possible_ins(first, third, ideal):
    return (third - first == (ideal - 1))
 
def is_possible_del(first, third, ideal):
    return (third - first == (ideal + 1))

INS = "INS"
DEL = "DEL"

def find_compare_positions(possible_first_pos, possible_third_pos):
    ideal_pos_diff = int(L*2/3)
    compare_positions = []
    for i in range(len(possible_first_pos)):
        for j in range(len(possible_third_pos)):
            if possible_third_pos[j] - possible_first_pos[i] == ideal_pos_diff:
                return []
            elif is_possible_ins(possible_first_pos[i], possible_third_pos[j], ideal_pos_diff):
                return [(possible_first_pos[i] + int(L/3), INS)]
            elif is_possible_del(possible_first_pos[i], possible_third_pos[j], ideal_pos_diff):
                return [(possible_first_pos[i] + int(L/3), DEL)]
    return compare_positions

def compare_insertion(ref,read_third,start_pos):
    diff = []
    num_compares = min(len(ref), len(read_third))
    for i in range(num_compares):
        if ref[i] != read_third[i]:
            diff.append([read_third[i], start_pos + i])
            break
    return diff

def compare_deletion(ref,read_third,start_pos):
    diff = []
    num_compares = min(len(ref), len(read_third))
    for i in range(num_compares):
        if ref[i] != read_third[i]:
            diff.append([ref[i], start_pos + i])
            break
    return diff

def evaluate_positions(compare_positions, thirds, ref):
    ins = []
    dels = []
    if len(compare_positions) > 0:
        for position_tuple in compare_positions:
            genome_to_compare = ref[position_tuple[0]:(position_tuple[0] + int(L/3))]
            if position_tuple[1] == INS:
                insertion = compare_insertion(genome_to_compare, thirds[1], position_tuple[0])
                if len(insertion) > 0: 
                    ins += (insertion)
            else:
                genome_to_compare = ref[position_tuple[0]:(position_tuple[0] + int(L/3))]
                deletion = compare_deletion(genome_to_compare, thirds[1], position_tuple[0])
                if len(deletion) > 0: 
                    dels += deletion

    return ins, dels

def find_insertions_deletions_in_read(read, lookup, ref):
    thirds = split_into_3(read)
    compare_positions = []
    if thirds[0] in lookup and thirds[2] in lookup:
        compare_positions = find_compare_positions(lookup[thirds[0]], lookup[thirds[2]])

    return evaluate_positions(compare_positions, thirds, ref)

def find_ins_dels(reads, lookup, ref):
    possible_ins = []
    possible_dels = []
    for read in reads:
        insertions, deletions = find_insertions_deletions_in_read(read, lookup, ref)
        if len(insertions) > 0:
            possible_ins += insertions
        if len(deletions) > 0:
            possible_dels += deletions

    insertions_to_count = count_occurences_possible_snps(possible_ins)
    deletions_to_count = count_occurences_possible_snps(possible_dels)

    ins = choose_majority(insertions_to_count, 10)
    dels = choose_majority(deletions_to_count, 10)
    return ins, dels

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
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

def write_lookup(lookup):
    file = open("lookup.txt", 'w')
    for key in lookup:
        file.write(key+"\n")
        [file.write(str(val) + " ") for val in lookup[key]]
        file.write("\n")

def read_lookup():
    file = open("lookup.txt", 'r')
    lookup = {}
    for key in file:
        lookup[key[:-1]] = [int(val) for val in next(file).split()]
    return lookup

def write_reads(reads):
    file = open("parsed_reads.txt", 'w')
    for read in reads:
        file.write(read + "\n")

def read_reads():
    file = open("parsed_reads.txt", 'r')
    reads = []
    for read in file:
        reads.append(read[:-1])
    return reads

def write_snps(snps):
    file = open("snps.txt", 'w')
    for x in snps:
        file.write(','.join([str(u) for u in x]) + '\n')

def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    lookup = create_subsequence_lookup(reference)
    # write_lookup(lookup)
    # lookup = read_lookup()
    reads = convert_pairs_to_reads(input_reads)
    reduced_size_reads = enumerate_reads(reads)
    # reads = read_reads()
    # write_reads(reduced_size_reads)
    snps = find_snps(reads, lookup, reference)
    # write_snps(snps)
    insertions, deletions = find_ins_dels(reduced_size_reads, lookup, reference)
    insertions.sort(key= lambda x:x[1])
    deletions.sort(key= lambda x:x[1])
    # insertions = []
    # deletions = []
    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

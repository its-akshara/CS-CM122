import sys
import argparse
import time
import zipfile

L = 30
MISMATCHES_ALLOWED = 4 # Number of mismatches allowed
CONSENSUS_MAJORITY = 20 # Number to get consensus that SNP located there

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

def create_subsequence_lookup(genome):
    subseq_to_index = {}
    for i in range(int(len(genome) - L/3 + 1)):
        seq = genome[i:i+L]
        if seq in subseq_to_index:
            subseq_to_index[seq].append(i)
        else:
            subseq_to_index[seq] = [i]
    return subseq_to_index

def split_into_3(read):
    return ([read[i:(i + L/3)] for i in range(0, len(read), int(L/3))])

def find_pos_differences(ref, read, start_pos):
    diff = []
    for i in range(len(read)):
        if read[i] != ref[start_pos + i]:
            diff.append([ref[start_pos + i], read[i] , start_pos+i])
    return diff

def calculate_ref_read_position(which_third, index):
    return index + which_third*L/3

def evaluates_indices(indices, read, ref):
    snps = []
    # for index in indices:
    #     ref_read = ref[index:(index+L)]
    #     diff = find_pos_differences(ref_read, read, index)
    #     if len(diff) < M:
    #         snps += diff
    return []

# returns list of [OriginalAllele,SNP,Position]
def find_possible_snp_in_read(read, lookup, ref):
    possible_snps = find_pos_differences(ref, read, 0)
    thirds = split_into_3(read)
    which_third = 0
    for third in thirds:
        if third in lookup:
            possible_indices = lookup[third]
            possible_snps += evaluates_indices(possible_indices, read, ref)
        which_third += 1
    # for i in range(1, len(ref) - len(read) + 1):
    #     diff = find_pos_differences(ref, read, i)
    #     if len(diff) < MISMATCHES_ALLOWED and len(diff) < len(possible_snps): # choose best case match
    #         possible_snps = diff
    return possible_snps

def count_occurences_possible_snps(snps):
    snp_pos_to_count = {}
    for snp_tuple in snps:
        if snp_tuple in snp_pos_to_count:
            snp_pos_to_count[snp_tuple] += 1
        else:
            snp_pos_to_count[snp_tuple] = 1
    return snp_pos_to_count

def choose_majority_snps(snp_possibilities_to_count):
    snps = []
    for snp_tuple in snp_possibilities_to_count:
        if snp_possibilities_to_count[snp_tuple] > CONSENSUS_MAJORITY:
            snps.append(snp_tuple)
    return snps

def find_snps(reads, lookup, ref):
    possible_snps = []
    snps = []
    for read in reads:
        possible_snps += find_possible_snp_in_read(read, lookup, ref)

    snp_possibilities_to_count = count_occurences_possible_snps(possible_snps)

    snps = choose_majority_snps(snp_possibilities_to_count)

    return snps

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
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
    reads = [read for read_pair in input_reads for read in read_pair]
    snps = find_snps(reads, lookup, reference)

    # snps = [['A', 'G', 3425]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

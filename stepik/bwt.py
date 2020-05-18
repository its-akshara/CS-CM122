L = 30
MISMATCHES_ALLOWED = 1 # Number of mismatches allowed
CONSENSUS_MAJORITY = 10 # Number to get consensus that SNP located there

def create_subsequence_lookup(genome, L = L):
    subseq_to_index = {}

    for i in range(int(len(genome) - L/3 + 1)):
        seq = genome[i:int(i+L/3)]
        if seq in subseq_to_index:
            subseq_to_index[seq].append(i)
        else:
            subseq_to_index[seq] = [i]

    return subseq_to_index

def split_into_3(read):
    L = len(read)
    return [read[:int(L/3)], read[int(L/3):int(L*2/3)], read[int(L*2/3):]]

def ref_start_pos(index, which_third, read):
    if len(read) % 3 == 0:
        result = int(index + ( len(read) * which_third)/3)
    else:
        result = int(index + ( len(read) * which_third)/3)
    print("Start at {}, given index {}".format(result, index))
    print("Read: " + read)
    return result

def find_pos_differences(ref, read, start_pos):
    diff = []
    for i in range(len(read)):
        if read[i] != ref[i]:
            diff.append(start_pos+i)

    return diff

def evaluates_indices(indices, read, ref, which_third):
    snps = []

    for i in range(len(indices)):
        ref_start = ref_start_pos(indices[i], which_third, read)
        ref_subseq = ref[ref_start:(ref_start+len(read))]
        # print("Compared ref: " + ref_subseq)
        # print("Read: " + read)
        if (len(ref_subseq)) == len(read):
            diff = find_pos_differences(ref_subseq, read, ref_start) 

            if len(diff) <= MISMATCHES_ALLOWED:
                snps += diff

    return snps

# returns list of [OriginalAllele,SNP,Position]
def find_possible_snp_in_read(read, lookup, ref):
    possible_snps = []
    thirds = split_into_3(read)
    which_third = 0
    print(thirds)
    
    for third in thirds:
        if third in lookup:
            possible_indices = lookup[third]
            possible_snps += (evaluates_indices(possible_indices, read, ref, which_third))
        which_third += 1
    
    return possible_snps

def read_lines(path):
    file = open(path, 'r')
    lines = []
    for line in file:
        if line[-1] == '\n':
            line = line[:-1]
        lines.append(line)
    return lines


def bwt_mx(text):
    bwt_mx = []
    rotated = text
    for i in range(len(text)):
        rotated = rotated[-1] + rotated[:-1]
        bwt_mx.append(rotated)
    bwt_mx.sort()
    return bwt_mx

def bwt(text):
    mx = bwt_mx(text)
    bwt = ""
    for rotated in mx:
        bwt += rotated[-1]
    return bwt 

def inv_bwt(bwt):
    inverted_mx =  ["" for i in range(len(bwt))]
    for i in range(len(bwt)):   
        for j in range(len(bwt)):
            inverted_mx[j] = bwt[j] + inverted_mx[j]
        inverted_mx.sort()
    
    for i in range(len(bwt)):
        if inverted_mx[i][-1] == "$":
            return inverted_mx[i]
    return inverted_mx[0]

def range_of_symbol(last_col, symbol, top, bottom):
    for i in range(top, bottom + 1):
        if last_col[i] == symbol:
            first = i
            while i <= bottom and last_col[i] == symbol:
                last = i
                i += 1
            return first, last
    return -1, -1

def create_last_to_first(bwt):
    first_col = sorted(bwt)
    last_to_first = {}
    number_seen = {}

    for i in range(len(bwt)):
        if bwt[i] in number_seen:
            number_seen[bwt[i]] += 1
        else:
            number_seen[bwt[i]] = 1
        
        last_to_first[i] = first_col.index(bwt[i]) + number_seen[bwt[i]] - 1

    return last_to_first

def bwt_matching(last_col, pattern, last_to_first):
    top = 0
    bottom = len(last_col) - 1

    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            top_index, bottom_index = range_of_symbol(last_col, symbol, top, bottom)
            last_short = last_col[top : (bottom + 1)]
            if top_index != -1:
                top_index  = last_short.index(symbol) + top
                bottom_index = len(last_short) - last_short[::-1].index(symbol) + top - 1
                top = last_to_first[top_index]
                bottom = last_to_first[bottom_index]
            else:
                return 0
        else:
            return bottom - top + 1
                
def read_patterns(path):
    file = open(path, 'r')
    patterns = file.read().split()
    return patterns

def write_positions(positions, path):
    file = open(path, 'w')
    output = ""
    for position in positions:
        output += str(position) + " "
    file.write(output[:-1])

def read_text(path):
    file = open(path, 'r')
    return file.read()

def find_exact_matches(patterns, lookup):
    positions = []
    for pattern in patterns:
        if pattern in lookup:
            positions += lookup[pattern]
    return positions

def find_matches_with_faults(patterns, ref):
    positions = []
    L = int(len(patterns[0]))
    lookup = create_subsequence_lookup(ref, L)
    for pattern in patterns:
        if L != int(len(pattern)):
            L = int(len(pattern))
            lookup = create_subsequence_lookup(ref, L)
        positions += find_possible_snp_in_read(pattern, lookup, ref)
    return positions


# text = read_text("input.txt")
# last_to_first = create_last_to_first(text)
# bwt = ""

# patterns = read_patterns("dataset_317416_4.txt")
# positions = []
# for pattern in patterns:
#     positions.append(bwt_matching(text, pattern, last_to_first))
patterns = read_patterns("dataset_317417_10.txt")

ref = "ACATGCTACTTT"
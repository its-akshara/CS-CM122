def read_lines(path):
    file = open(path, 'r')
    lines = []
    for line in file:
        if line[-1] == '\n':
            line = line[:-1]
        lines.append(line)
    return lines

def add_pattern(trie, pattern):
    current_node = 0
    max_num = max(trie.keys())
    for i in range(len(pattern)):
        symbol = pattern[i]
        if symbol in trie[current_node]:
            current_node = trie[current_node][symbol]
        else:
            max_num += 1
            trie[max_num] = {}
            trie[current_node][symbol] = max_num
            current_node = max_num
    return trie

def construct_trie(patterns):
    trie = {} # map node -> (letter -> next_node)
    trie[0] = {}
    for pattern in patterns:
        trie = add_pattern(trie, pattern)
    return trie

def is_leaf(node, trie):
    return len(trie[node]) == 0

def prefix_tree_matching(trie, text):
    index = 0
    symbol = text[index]
    node_num = 0

    while index <= len(text):
        if index < len(text):
            symbol = text[index]
        if is_leaf(node_num, trie):
            return True
        elif node_num in trie and symbol in trie[node_num]:
            index += 1
            node_num = trie[node_num][symbol]
        else:
            return False
    return False

def find_pattern_matches_in_trie(trie, text):
    matches = []
    i = 0
    while len(text) > 0:
        if prefix_tree_matching(trie, text):
            matches.append(i)
        text = text[1:]
        i += 1
    return matches

def write_trie(trie, path):
    file = open(path, 'w')
    for node in trie:
        for edge in trie[node]:
            file.write("{}->{}:{}\n".format(node, trie[node][edge], edge))

def write_positions(positions, path):
    file = open(path, 'w')
    output = ""
    for position in positions:
        output += str(position) + " "
    file.write(output[:-1])

def construct_noncompressed_suffix_trie(text):
    trie = {} # map node -> (letter -> next_node)
    trie[0] = {}
    max_num = 0
    leaf_to_label = {}
    edge_to_position = {}
    for i in range(len(text)):
        current_node = 0
        for j in range(i, len(text)):
            current_symbol = text[j]
            if current_symbol in trie[current_node]:
                current_node = trie[current_node][current_symbol]
            else:
                max_num += 1
                trie[max_num] = {}
                trie[current_node][current_symbol] = max_num
                edge_to_position[(current_node, current_symbol, max_num)] = j
                current_node = max_num
        if is_leaf(current_node, trie):
            leaf_to_label[current_node] = i
    return trie, leaf_to_label, edge_to_position

def compress_suffix_trie(trie, leaf_to_label, edge_to_position):
    return trie

def construct_suffix_trie(text):
    trie, leaf_to_label, edge_to_position = construct_noncompressed_suffix_trie(text)
    return compress_suffix_trie(trie, leaf_to_label, edge_to_position)

# def MaximalNonBranchingPaths(Graph)
#         Paths ← empty list
#         for each node v in Graph
#             if v is not a 1-in-1-out node
#                 if out(v) > 0
#                     for each outgoing edge (v, w) from v
#                         NonBranchingPath ← the path consisting of single edge (v, w)
#                         while w is a 1-in-1-out node
#                             extend NonBranchingPath by the edge (w, u) 
#                             w ← u
#                         add NonBranchingPath to the set Paths
#         for each isolated cycle Cycle in Graph
#             add Cycle to Paths
#         return Paths

def calculate_node_degrees(trie):
    in_degree = {}
    out_degree = {}

    for node in trie:
        outgoing_edges = trie[node].keys()
        out_degree[node] = len(outgoing_edges)
        for edge in outgoing_edges:
            if trie[node][edge] in in_degree:
                in_degree[trie[node][edge]] += 1
            else:
                in_degree[trie[node][edge]] = 1

    return in_degree, out_degree

def is_one_to_one(node, in_degree, out_degree):
    return node in in_degree and in_degree[node] == 1 and node in out_degree and out_degree[node] == 1

def maximal_nonbranching_paths(trie):
    paths = []
    in_degree, out_degree = calculate_node_degrees(trie)

    for node in trie:
        if is_one_to_one(node, in_degree, out_degree) == False and out_degree[node] > 0:
            for outgoing_edge in trie[node]:
                print("not one to one and not leaf")
                nonbranching_path = [(node, outgoing_edge, trie[node][outgoing_edge])]
                current_node = trie[node][outgoing_edge]
                outgoing = list(trie[current_node].keys())[0]
                while is_leaf(current_node, trie) == False and is_one_to_one(current_node, in_degree, out_degree):
                    outgoing = list(trie[current_node].keys())[0]
                    nonbranching_path.append((current_node, outgoing, trie[current_node][outgoing])) 
                    current_node = trie[current_node][outgoing]
                paths.append(nonbranching_path)
    return paths



# trie = construct_trie(read_lines("dataset_317406_8.txt"))
# write_positions(find_pattern_matches_in_trie(trie, "CTTTGATACTGGGCCGCCTGCGGTTGCGCAGGTCGGTAGGTGACCCTATTTGTCCTTAAGCCCCATGGCATTAGTCGCGGATCTAAAGGGTTGGGGGGCCGTGACACTCGTCGCTATTTGGGCACAATCCTCAACTCGAGGTCGTTCTACTTCTGTGGTGTCCCATCAGGCATTAGCATATGTGAACTTACCACTGTATTCCTCTTACGGGGTGTCAGTCTAGTCCTAGCCTCGACTGGTTATTTATACGAGCATTTTAGGGCGGAGGTTACGAACCGCACGATGTCTGACCTATCATCCGGAGCTGAGCACCGACGCTAGGCGAAGCACTATGCCTGCGTGACGTGCCGGCCGCAGCGTCCTTGGCTCCGTCGGACTCTATTGGCTGGACTCCTTCATCATATGTTCAGTGCGTCTGTTTTGTAGCAAAGACCCCGTACTAAATGCGCTACCGATCTAACGCTAAAGAACGAAGTGCATTGAGGGATCGGGGGATAGGCCCGTACATTTGACTCACAGCGAGCCGACCTATTCAATCCCGGTGCCTCCAAGAGTAACTGCAATTAGCCGATAGCTTCAAGACACGTTTACCTCGCGAACTATTATTACCGATTAATCACCAATCCCGCAACCGTCGTCGTCACGACTCGTACGAAGCGAATATGCGCGGAGCTCCACCTGTTGAATATGCGAAACGAATGAAACATCAGTCTGGTCGTTCGAGCGCGGCGGCCGGGCAATGAGAATCGCGTGCTGATCACTAGCAGTTGGTAGGGCCACCGTTTATTAACTACGGCAAATGGCGTCCGAACGGCCACGTGTCTTGCATCTGCATACAGGGATACGTCTGGTTGGCGGAGAATTTCATGAGACAATATCAAAACCTAGTGTCAAGAAGGGCGGGACGTAGAAAGGTAGCGTCCGAAATCCTGGAGCTGCGCGATTACCGATTACCGCATCTCAATGTCCGTAACCACGAGAATGTGCCCGACAGGGGTCAACCTTGCGTTAAACCTTTGTGATTATTCGATCATCTGCACTTCGTATAGGGTGCCTAGCTAGTGCAGGAAATTTTGCACTCTCGAAACCGGATAAGACCCTTAACCAGGTGCTGGTTGAGTGATCACGCAGCTCCGGCTGAGATTCAAAGCGCACACAGGTCGTAGGATAGGGATTAGGTATCGAGAATAACAAGGAAGATACTATGGGAAGACCTCAGACTTATGAGTAAGGGAGCAAAAAGCAGGCGGTCTGCGCGCCCCAGTGCAACCAGACTTGTCTTTCGTATTCTGAAACTAATGACTCGACGGACCTGAAGGGCGTGGATCAGCCATAGTAATTAAGAAAGACACCACTACAAAAAGCGAATTGTTATGTGCGGCCGACGACATGCCCATAAATGATTGATGCTTAACTTAAGATCATATGCGAGAAACAGTAGTTGAACGGCTTCCTCTGGAACTGAGGGGTCTTATAGGCGGCCGACGTAGGGCACAAACAAGAATCAGTTGCTTCTCCATTTCATAAGCATGAATTTACCGAGTGCCATACACGAAACACCACGTGACCAAAGTTGTTCTCCATTCTCCATTTCTAGCTTTTCCCACGGGACTTTCCGACCGGGAGTGGCATCATTCCTTCGGAAATTGAAAGAAATGCTACGTGGTCCTTTCCAATCCAGTGGGAGGAGTGCTATTATAGCCTGAAAGTAAGGTACTTACACGCGGGTAAGCTTATCCGGGTGTTAGTGTTCGTGAAGTTGGGCCATAGGGATTCCGTTGTGTCGTCAGCATATTTATGTGAGTCTTCACTAGCCGCTGGATTTTCATTACGCGGGGTGGATACCATGGAGAATGCCCGCGCGTCTGCGGTGACAGGTAAAGATGAGATCAACTGCCAACTGCCATCTGCTCGCATTTGCTTGTCTCTAGCGCTAGACATGGTCATCCGATTATCCCCGAAGGAACCGAGATTAACAGGACTCGACTGAATCGCTTACCGCAGCCGGCAGAGCGGCACAGTCTAAGACCCTTTAAATCTTCGCCGGGCCGGAGAAACGGGGCAACTGCCAGGCTCTGCCTCATGACACGTCAAATAATAGCAAGAGTAGGTGATGATGCATGTGAGAATCATGAGAGGTGTAGCACTGCTCGATACCTGATCTATAGCGGGCATGTGGCTGGGCTCACTCGACTTATGACATACCCTTTCAAAACCCACGAGCAAATCTTTCAGGGGTTCGTCTAGGTAACGCAACAGGGTGAATCAGATGGAACAGTATTAACACCGCCATATATAGCCCTTAAGCCTTAAGCCTTTCTCCATTTTAACCGTGTATCCCGATGTTGAGTAATGTACTTCACTGCGGTGGAACGTAACAGCGACAATTCCATCGTTCGTCGAGAGCCGCCTCGAATCTATCTCCCTTATGGATTCGTCTGAGCCTAACGACCATCTTAGGTAGTTGTCATTGACGGCTTAGCTCCGGGGACGATTTTCTCCTGGGGTTCGCGGGAGAAGTTAGGCAGGTACAAACTTCGTTCTGCACATCCAGCCATCTCGTGGTATGTGGTATGTGACTGAATTATTGTACAAGCCTGCGTTATAAGTGTAGGGATAGCCCGAGTAGGGACACCACGGGTACCATCTCAGTTCTAATTAGGGGTCTCCCGACCAACGTAGAAGTCAAGAAGGAGCCACGGCGAAAATCTGCTGGGCCTGTACGGTGGAAGCAGATGGTGTCTCGCAAGGACGGACGCTCCTTCTCCAAACATAGTATTAAGTAGGTGCAGATAATGAATCTAGAAGGGGTCCCTAGTGTCCCCACTCTCTACGCCATTCTAAAGCTACTGGCATCGCTTTACAGCCGGTTTTTCAAACCTTAAGCCTTAAGCCTGCGACCTACGCCTGAAGTTTTCACTCGTGCCTTAGGTTGGTGACTTCGGAGGATAACATTGCTTTGTACCCGTGTTTCTGTCTTCCCCTTCGGTACAGGTGGCGAGGCAGGTCTGAAAGCCGCGACCCTTTTACAAAGCAACTTAATATTGGATGGGTGTAACAACGCCCAGCTTAGTCCGCCCAGCTAGGTGTCCAAAAAAAGCGACAATCGTCACGAATAGGAAGCAGCCGGATTCTCCGTAAAGTGCAATTCTCACCTCCATTGGTACCGATAATAAATATGAATTCAGTACGGGCCCTGAGTTGGCCTTGGTACTGGCCATCTCTATCCTGTTAAAGTGCGAAGGACTATGCGTTCGCCGGAAGGGGTGGCACGAACCCTAGGGCAGTGGAAGCTTAGTCGGTAAGCTGGGCAGAGATGGAACGGTATGCACTGCTAACACGGCTATTGTGCTTGTGGTACAGCTTCACGAGGCCTGTAATCTGTGGACGGACGTAGATGTACTTTGCCGCTTTTGGTCGCCCGACGGACAAAATAATGCCAATCTCATATTTTTGGCAGCATACCCATATCGTTACTCCGACCCCGCGGTGCTCCTGAAAGTCTCTTAGCCTATTAATAAAAGTTAGCCTAGAAAATGGTAACTAGTATCTGAGCTTGCTAGCGATACCAAGTCGCCTCCGGAGGCAGAGATACCAAATGAATTACGTGGTAGCAGGGCGCATCATAACTCCATCAAGCGCCCCGCACGGTCCTTAATAACTCGCCGTGTTAGAATTACCATTGGCTGTAGTACGACAATCCCTTCGCCATGTATTAGGAATGACGATTCCCGGGAACGCATAATAGGAAAGGTCAAAAGTCTGTCGGTACCAACCTGCCGCCAAGTTTTTCGTTTCCATTGCTGCATAATTGTGGAATATGCATCTAGGATTCTAGACTAGATTTAGTGGGGAACGGCAGGGCGCAGGGCGCAAAGTCCGGCCGCTGCTGCTTAGCAGGAGCACACTAACCCGGTGTCTGCCGTACCTAGTCAACTTCCTCGCGGTCCTCTTAGCTCATTGCTGTCAGCGCGCACCGATTGCGAGTTGATCAGAACTTGTCCTTGCGCAGTTCAAAACCCGAGCCTGATCTAGGATACGCGTGCCCAGATAGGGCTTCTTTCCATTCTTCATGGGACTACGCTATCTTTCGCCATACGATTGCCGTGGCTTTACCGTCTGGGGCTAAAAGACTGAGCGTACTCAGGACACTCTAGTTTGTACTTTCCGGGTATAATCCTCCGGCGTTTAGGCTTGTTTAGGTACGGGTATGAACACCAAAACTGATCTCAGTTCTCAGTTCGAAATCGACCTGACGTATCTTGCGGCCCTTATGCAGATCGGCAGCCCAAGTTCCGCTCCAAGGAACATTGCACCGTCATGTTCATACCTTCATAATACGTGTTGACCTCAGCTTTAAGAGCCATGTGATTCATGAACTCCTACCCGATGCCTAGCAGGTTGGCATAATAACACGGGTTAGTTATGGTATGTGGTATGTGAACACAAGTCTACAACTTGATTGCTGTACTGTTTGGCACATGAGGATATATGCGCCCGTGGAGACGAACCGGAAAATTCGCATGGTATGTGGTATGTGTGGGGGGACGGCTGTTAACAGCGAAAACGAAAGATTTTAGGAAAGAGTTCTAGATGTAACATTCGGCGTGCCGCTCACCCGATGTTAAACCTGCAAGCCGTCAGCCACATGTAGGATTACTCTGACGGGTGGAGTGCGTGATCTCCAGTGCCGGGTAAAATTGCGCAACTGCCAACTGCCAGGGGATGATCTTAGCGCGGAAAGCTCAAGCAGGTTCACTAGTATTCCTTACGTGCCCCTCGAGTCAGAAACCCCGTACTTTAAGGCAGCGGTCGATTGGCTGGCCTCGCTTTGTCTGTGAGCAGGGCGCAAAGTAAATAGGTGCTCCGACCCTTTGTATGTGCCACCCGCCTCATTGAAGGGTGATGGTCTAACTTGCTTATCCAGTAAACCGATTATGATTGATAATATGAGCAGCCGCCGCCGCTTATTGTTGACGTTCAATCAGTCCCTCGCTGTCGTGGGGGAGTACACTATGGGGTATTCACCAGAAAAATTGTTATATCGTGATGAGCTGAGATTAACCTCTTATGATTGCAAGCCGTAAAGTTCAACTTGTTGGTATGTGGTATGTGTGCGCGTAGGCCCTTGGGGTATTGTCTAGAAACCTGATGCGGTGGGAACACATCCGACGTGGCAGCGTTCGACTGTCATAGAGTAAGGCACACTGTACGCCGACAGACACCATCGCTCTCACGTCTGGGTCGCAGGGCGCGGATCTTTGCACTATAGTTGTAAATACCGTATACTCAACGTCGTGATCAATGGGTCCAACGCCTATTTACCACTATGAAGTGTATTGTTGAAGAAGCGAATGACGCTTACTTCGTTTCCTTTCTATACTATAGGTATAGTCAGATTAGATCCGTACGGTACGTCAAGGTGCGTCACAGGCGCCCTAAGTACAAGCTGATAGTGTTTGGTTAATACTTTGGGCACGTAGCCAAAGGAGGGTACAATCAGTGTGGACTCTCCGAACTTAATAACGTTATAGGTGGTCTCTCGACCTTTCCTTTCACTGAGCTACACCTGGTCGTCGGGTAAGAAGCGTAATATATATCGCTACTTACACACGGCTTCGGCGCCGTACTCAGACGGAAGTCAGATCAGAAGACTTATTCGTGTGTTTTTCGAGTGCGTAGGCGCCCGCTGGCGATAAAATTAATCTCATACAAGACTCGATTGCCGTTAGGGCGCTAAGACTCGGGGTCGGATTTACTCAATACGGGGCTATGAGGTGAAAGAGTAGGTGGGGTCTCAGTTCTCAGTTCTGACTGCGGCACAAGTTTCACCTCCTCGGTCCCGTTTGTCACCTCCCTTGCTCAGCGCAATAACGTTTCCAGAATCTACAGTCTAACGGAAAGCGCTCTACAATCCCTGTCGATTACCGCTCCAGGCTCAGCTCACGTTCTGGAGTGCATAGGCCTATGTCGTTAAAGACAGAATAAAAGTTTCGTGACGAGGGGTTCATGACTCATGAAAACGAAGTCTAGGGGTTCTCAGCGTGTTACCCGAGAAGATACCCAGGCCAGTAATAGCCTCCACTTTTCGCCCGCTGATATAACCCGTCTTGCTCCAACTGCCAGTTCACGATCTGCCGCGAGGTATGAATAACCGTCGTGGGCACTCAATCGTGTGAGGTATGAGCTCAGTTAACCCCTTCCTCATGCCTAACACGTCCAAAGAAAAGATGGGACCCCGACCGTCACCCCGAGGCCGTCACGCGCCCCCATTCCCGGGGAAATTAACATGAAGTCTCCGCTACTTTAGTAGCTTGTACCAAGGCCGTAGCCCGATTACCGTCCGCGCAGGGCGCAGGGCGCAAAATGATTTACGCGAGAGGCTTGGCATCCGCACCGAAGCGGTCTAACTTTAGCCGGGCCCAGCATATTGGCAAGTCACTTTCGTCGGAGAAGCAGGACCCCGTACCGCCAGTTTGTGACACATGATTCGGACGCCATTGCGTCATGCGCTGACTCCAATAACGTCATCTAACACAATTTGGATATTGGCTGCGCGAACATCAAGACCCCGTACCGTGTCGCTTCTCCTGTTACATCTTCCCGATTACCGATTACCGCCGTAAAATTAATAGCAAATACGAGGTGCGAGGACACAATCTCAACATAATGTCGTGCGTAGCCACTAATTCTCCCCCTCCCCAACGGTGGGGGGCCACGGTTCACAGCCACGATTACCGGGGCTTACCCCGCGCTCACCGGTCCCGAATTAAGTCAATGCGACATAAGCGTTGGAAGCACGGCTATGCTAGCTGCTTAGTGGGCCAAAGCAAACACGACGTGAGGCCGCATTCGGGAAGACACGGCCATCTGGATAGAAGCCAGGGGAAACGTCTGTTTCGAAGTGGCTCTCGGCTTCCAGTCAGTCATGATCTGTCCTCGCCGTGAGATCAATCCACCTCGGTGAACTGGACCTCTATATATCTTCTCCATTCTCCATTCGCGGCTCCCTGCCGCGGACTGAAGCTAGTCTCACGATCCTAGGTCCTGGATTCGTGCGCAACTCTACCTCGAAAGATCGGAGGCCGAGTAGGCCTAATTCATTAGTAAGTTGAAACCGAATGAGGAACATGTTGCGTTCCTGTGGTATCAATCGGGAGCCCAGTATTTCCCACTCGGAAAGCGGGGTGTCGACGTGGGTTTTACGTACAAGACCCAGGCTACCGCAATGTAGCCTAAGGTACATGCGGTTTTAGTGAATATTTACATACGTTCCCCAGTGAGAGCCTGGACGAATGCCACTACGATCCTAGTCTTTCCTCAACGAGCTGGGTACTAACAATCCTACGTCATCTTCATACCGATGCTTCACACAGAGAACCAGGACTTGCGCAAGGGATGGGACGACGGCAGCGTAAGACAACCATTCTTGATTGCACGAGGACCGAGTGACTGCCGGGCGTTCCAAACGACGATGGCGTGTCATGTCGGCTCCCAGCCCCTCCTTTCTTAAGTCTTTAAAGAACGATGTCCTTCGACTAGTTGAATGGGCTTCTCGGACTCATCCGCTAGAACGTACAGGAGGATAAAGGACAGATCGAGTCAGTTACCCGTCTAGAGTACCCCTCTAATCCCTCTCTGACCATGTACGTACGCTCACATCAATGTTGCCGTCGGGTACAAGAGCCTCCCTTAAGCCTAATGTCGACAGGAGGCTACGCGTGCGGCGGTGAACGTTCAAGTCGGTCTAGCCCATCAAGTACACGTTTACAGCGCCCGGTTAGCGACGACAGGTTCTTTCCTATCCTAATGGAAGTTTATGGAATAAATTCGGGCTGCGAGGCCACTACCCGTATTAATACTGTGGTCTTTGCCACTTTTCTTTCGAATCGATCTCACGTAAGAACTGGGTCCTGAGTCCGAGGCGTAGGGGCGCGGGCCCTATAGCTACCCGCTGTCGAATTACCTAGCTGTCCAGCACAGTCCACAACAGTTCTCGTGGTGGCATGTTAATACCATGGTGATTGCGAGGTACAGAACAGGGGCTTATGCAGCCTTAGCAGTTATTTCTCAGTTCTCAGTTCTTGTCTAGACATCGTTTAGATCTGTGATCGTTCAGTAAGAATTTGATGGCGAACCTCTTCCCCCGTAGAAGCCACAAGTGGTCGGGATAGGGGCGGAGGAGCTCCGCAGGCGTTGAGACCAGCTCGAATGGTCCGAAATAGTCCTTTGGTCCCCGGGGGGGGCGATAAGTCACTTGGTGTCGCTCGCGTGACCCATCTATTTGCCTCAAAGAACGGTCACCATCCTGCTCCTTGCAACTATTATAAATAGGTCCAGTGTAGATCAAGTAGGGCGTACGTCTGCGAACTAATTCGTTAGCA"), "output2.txt")
trie = construct_suffix_trie("paaa")
print(trie)
print(maximal_nonbranching_paths(trie))
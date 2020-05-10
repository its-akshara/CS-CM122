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

    while True:
        if is_leaf(node_num, trie):
            return True
        elif node_num in trie and symbol in trie[node_num]:
            index += 1
            symbol = text[index]
            node_num = trie[node_num][symbol]
        else:
            return False
    return False

def find_pattern_matches_in_trie(trie, text):
    matches = []
    i = 0
    while len(text) > 0:
        if prefix_tree_matching(text, trie):
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
        output += str(position)
    file.write(output)


print(find_pattern_matches_in_trie(construct_trie(read_lines("sample.txt")), "AATCGGGTTCAATCGGGGT"))

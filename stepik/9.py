def read_lines(path):
    file = open(path, 'r')
    lines = []
    for line in file:
        if line[-1] == '\n':
            line = line[:-1]
        lines.append(line)
    return lines

print(read_lines("sample.txt"))
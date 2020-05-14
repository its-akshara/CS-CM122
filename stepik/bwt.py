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

                



# print(bwt("TTCGGCGCCGGGGATAGAGTGTAGGTTACCCAATGCCCGACTCATACAGCGTGTGCGACTCTCCAACCGAGTCGTTAGTTTTTGGTCCCTGTCCAAGTTGCGAAGAAACAGTCCCCCTTAGCTATGCTACTAGCCCTGAGCCCCGGTATCGAAGAGATTAAGCGGAGAATCAAGCTCGTGGAACTATGACCTAACACCTCAGTCGTGATAGGGTCGAAGAAATAAACGCCCAGGAATTTCACTTAGTGATGCGCCACGGGGCAATTCCGCTGCATCGAAAAGTGTGTCTCACGCTTACCTCCACTAATAAAGATCTATCTGGTCATATCTAGATCCCACAGAAGCCTAGAGTCGAGGTTGAGTCGGGAAAAGGGACTGCTGAGGCGGTTGCCGCTTAGTGGGCTCCCAAGCAGCGGGGAGGTATTAAACTGGTACCGGGTTAGTCGAGTATGGCCATGGAAGCTTTTGATATCCGGATTTTGCGGCCGGGGAGGTGTACCAAAGTACCAGAACTATAGCAACGAAAGCTCAATGTAATCGTACGGCACCGGATCGCGGGCTCCCCTAATTGTACCACCGCACCTCCCCTTCCAGTTATGGGACTGCGAAGAGGCTCCACGCCTCGAGGGTTTGTATTGCCTGCTAAAGCTTATTTTCTATTCCAAGCCAAATCGGGTTAAGAGCAGTGCCTGTACCCGATTGCTCCTCGTGTGCCCGGGAAAGCAGGTGGAAAGTAATCAGCCTACCTAACCAGTGTAATTGTTCGGACACAATACGCCGAACGCGCTTGTCTTCATGTTAGAGCCATGTTCTTGTCTAACAGTCGATCGAGTTCGTCATGCTGCGCATCGGGGTGAACTCTACCAAACTACCCGTGATAGGCCGAGTCAGTGCCGCCCGACATCT$"))
text = "abcb$"
bwt = bwt("abcb$")

patterns = ["b"]
for pattern in patterns:
    print(bwt_matching(bwt, pattern, create_last_to_first(bwt)))
print(text[:-1])
import numpy as np

# def valid_coord(coord, max_coord):
#     return coord[0] >= 0 and coord[0] <= max_coord[0] and coord[1] >= 0 and coord[1] <= max_coord[1]

def count_paths(current):
    x, y = current[0], current[1]
    if x == 1 or y == 1: 
        return 1 
    return count_paths((x-1, y)) + count_paths((x, y-1)) 

def dp_change(money, coins):
    min_num_coins = [0]
    for m in range(1, money + 1):
        min_num_coins.append(999999)
        for i in range(len(coins)):
            if m >= coins[i]:
                if (min_num_coins[m-coins[i]] + 1) < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m-coins[i]] + 1
    return min_num_coins[money]

def read_int_list(path):
    file = open(path, 'r')
    result = []
    for num in file.read().split(','):
        result.append(int(num))
    return result

def manhattan_tourist(n, m, down, right):
    s = np.zeros((n+1, m+1))
    s[0][0] = 0
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + down[i-1][0]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] + right[0][j-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])
    return s[n][m]

def read_manhattan(path):
    file = open(path, 'r')
    n, m = int(file.read().split())
    return n, m

read_manhattan("input.txt")
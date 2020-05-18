# def valid_coord(coord, max_coord):
#     return coord[0] >= 0 and coord[0] <= max_coord[0] and coord[1] >= 0 and coord[1] <= max_coord[1]

def count_paths(current):
    x, y = current[0], current[1]
    if x == 1 or y == 1: 
        return 1 
    return count_paths((x-1, y)) + count_paths((x, y-1)) 

print(count_paths((17,13)))
import os

depth = 1000
ofile = open(f'../test/data/gfa_files/worst_case({depth}).gfa', 'w')

tot_nodes = 4 * depth
tot_edges = 6 * depth - 2
ch = 26

dict = {}

for i in range(tot_nodes):
    s = '' if i > 0 else 'A'
    x = i
    while(x > 0):
        s += chr(65 + (x % ch))
        x //= ch
    dict[i] = s
    print('S', s, '*', sep = '\t', file = ofile)

for i in range(depth):
    print('L', dict[4 * i], '+', dict[4 * i + 2], '+', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i], '+', dict[4 * i + 3], '+', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i + 2], '+', dict[4 * i + 1], '+', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i + 3], '+', dict[4 * i + 1], '+', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i + 2], '-', dict[4 * i], '-', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i + 3], '-', dict[4 * i], '-', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i + 1], '-', dict[4 * i + 2], '-', '0M', sep = '\t', file = ofile)
    print('L', dict[4 * i + 1], '-', dict[4 * i + 3], '-', '0M', sep = '\t', file = ofile)

    if i > 0:
        print('L', dict[4 * (i - 1)], '+', dict[4 * i], '+', '0M', sep = '\t', file = ofile)
        print('L', dict[4 * i + 1], '+', dict[4 * (i - 1) + 1], '+', '0M', sep = '\t', file = ofile)
        print('L', dict[4 * i], '-', dict[4 * (i - 1)], '-', '0M', sep = '\t', file = ofile)
        print('L', dict[4 * (i - 1) + 1], '-', dict[4 * i + 1], '-', '0M', sep = '\t', file = ofile)

ofile.close()
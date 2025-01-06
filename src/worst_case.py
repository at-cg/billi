import os

depth = 10
ofile = open(f'../test/data/gfa_files/worst_case({depth})_h.gfa', 'w')

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
    print('S', s, '*', 'LN:i:495', 'ng:i:33', 'nc:i:39', 'c1:i:33', 'c2:i:0', f'pp:Z:{s}:ENSP00000496625.1', sep = '\t', file = ofile)

for i in range(depth):
    print('L', dict[4 * i], '+', dict[4 * i + 2], '+', '0M', 'ng:i:33', 'nc:i:33', 'ad:i:8451', 's1:i:2181', 's2:i:1325', sep = '\t', file = ofile)
    print('L', dict[4 * i], '+', dict[4 * i + 3], '+', '0M', 'ng:i:33', 'nc:i:33', 'ad:i:8451', 's1:i:2181', 's2:i:1325', sep = '\t', file = ofile)
    print('L', dict[4 * i + 2], '+', dict[4 * i + 1], '+', '0M', 'ng:i:33', 'nc:i:33', 'ad:i:8451', 's1:i:2181', 's2:i:1325', sep = '\t', file = ofile)
    print('L', dict[4 * i + 3], '+', dict[4 * i + 1], '+', '0M', 'ng:i:33', 'nc:i:33', 'ad:i:8451', 's1:i:2181', 's2:i:1325', sep = '\t', file = ofile)

    if i > 0:
        print('L', dict[4 * (i - 1)], '+', dict[4 * i], '+', '0M', 'ng:i:33', 'nc:i:33', 'ad:i:8451', 's1:i:2181', 's2:i:1325', sep = '\t', file = ofile)
        print('L', dict[4 * i + 1], '+', dict[4 * (i - 1) + 1], '+', '0M', 'ng:i:33', 'nc:i:33', 'ad:i:8451', 's1:i:2181', 's2:i:1325', sep = '\t', file = ofile)

ofile.close()
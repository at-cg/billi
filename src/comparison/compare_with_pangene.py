import os

data = ['MHC', 'C4-90', 'human100', 'human100p10', 'Mtb152m-p0a1', 'Mtb152m-p1a2', 'Mtb152p']

for i in range(len(data)):
    panbubble_path = f'../../panbubble_data/results/small/panbubble/{data[i]}/panbubble.txt'
    pangene_path = f'../../panbubble_data/results/benchmark/pangene/{data[i]}.txt'
    out_path = f'../../panbubble_data/results/small/panbubble/{data[i]}/compare_with_pangene.txt'

    def get_string(g1, g2, s1, s2):
        if(g1 > g2):
            g1, g2 = g2, g1
            if s1 == '>' or s1 == '+':
                s1 = '-'
            else:
                s1 = '+'
            if s2 == '>' or s2 == '+':
                s2 = '-'
            else:
                s2 = '+'
            s1, s2 = s2, s1
        else:
            if s1 == '>':
                s1 = '+'
            elif s1 == '<':
                s1 = '-'
            if s2 == '>':
                s2 = '+'
            elif s2 == '<':
                s2 = '-'
            
        s = g1 + ' ' + s1 + ' ' + g2 + ' ' + s2
        return s

    lb = []
    cnt = 0
    with open(panbubble_path, 'r') as file:
        for line in file:
            tokens = line.strip().split()
            g1 = tokens[0][1:]
            g2 = tokens[1][1:]
            s1 = tokens[0][0]
            s2 = tokens[1][0]
            lb.append(get_string(g1, g2, s1, s2))     
    lb.sort()

    lp = []
    with open(pangene_path, 'r') as file:
        for line in file:
            tokens = line.strip().split()
            if(tokens[0] == 'BB'):
                g1 = tokens[4][1:]
                g2 = tokens[5][1:]
                s1 = tokens[4][0]
                s2 = tokens[5][0]
                lp.append(get_string(g1, g2, s1, s2))
    lp.sort()

    # checking for duplicates
    for i in range(1, len(lb)):
        assert lb[i] != lb[i - 1]
    
    for i in range(1, len(lp)):
        assert lp[i] != lp[i - 1]

    with open(out_path, 'w') as ofile:
        print(f"pangene count : {len(lp)}", file = ofile)
        print("Bubbles marked by pangene only:", file = ofile)

        for i in range(len(lp)):
            if lp[i] not in lb:
                print(lp[i], file = ofile)

        print('====================================================', file = ofile)

        print(f"panbubble count : {len(lb)}", file = ofile)
        print("Bubbles marked by panbubble only:", file = ofile)

        for i in range(len(lb)):
            if lb[i] not in lp:
                print(lb[i], file = ofile)

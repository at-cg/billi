import os

# data = ['Ecoli']
data = ['human100', 'human100p10']
# data = ['MHC', 'C4-90', 'Mtb152m-p0a1', 'Mtb152m-p1a2', 'Mtb152p']

for i in range(len(data)):
    billi_path = f'/home/daanish/projects/Pangene/results/Billi/{data[i]}/summary.txt'
    out_path = f'/home/daanish/projects/Pangene/results/Billi/{data[i]}/compare_intra.txt'

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
    cnt = -1
    with open(billi_path, 'r') as file:
        for line in file:
            tokens = line.strip().split()
            if tokens[0] == 'Total':
                cnt += 1
                lb.append([])
                continue
            if tokens[0] == 'Tips' or tokens[0] == 'FP':
                continue
            g1 = tokens[0]
            g2 = tokens[2]
            s1 = tokens[1]
            s2 = tokens[3]
            lb[cnt].append(get_string(g1, g2, s1, s2))     
    lb[cnt].sort()

    # checking for duplicates
    for cnt in range(2):
        for i in range(1, len(lb)):
            assert lb[cnt][i] != lb[cnt][i - 1]

    with open(out_path, 'w') as ofile:
        for i in range(len(lb[0])):
            if lb[0][i] not in lb[1]:
                print(lb[0][i], file = ofile)
import os

data = 'MHC'
in_file=f"/scratch/projects/daanish/data/Bubbles/Data/{data}/{data}.gfa"
out_file=f"/scratch/projects/daanish/data/Bubbles/Data/{data}/{data}-modified.gfa"

with open(out_file, 'w') as ofile:
    with open(in_file, 'r') as ifile:
        for line in ifile:
            print(line.strip(), file = ofile)

            if(line[0] == 'L'):
                tokens = line.strip().split('\t')
                g1 = tokens[1]
                g2 = tokens[3]
                
                if(g1 != g2):
                    s1 = tokens[2]
                    s2 = tokens[4]

                    if s1 == '+':
                        s1 = '-'
                    else:
                        s1 = '+'
                    
                    if s2 == '+':
                        s2 = '-'
                    else:
                        s2 = '+'
                    
                    s1, s2 = s2, s1
                    g1, g2 = g2, g1

                    tokens[1] = g1
                    tokens[3] = g2
                    tokens[2] = s1
                    tokens[4] = s2

                    t1 = ('\t').join(tokens)
                    print(t1, file = ofile)

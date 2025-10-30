in_path = ['/home/daanish/projects/panbubble_data/Data/CHM13/v2/chromosomes/chrY.gfa', '/home/daanish/projects/panbubble_data/Data/CHM13/v2/chromosomes/chrX.gfa']
out_path = ['/home/daanish/projects/panbubble_data/Data/CHM13/v2/chromosomes/chrY-modified.gfa', '/home/daanish/projects/panbubble_data/Data/CHM13/v2/chromosomes/chrX-modified.gfa']

for i in range(2):
    with open(out_path[i], 'w') as ofile:
        with open(in_path[i], 'r') as ifile:
            for line in ifile:
                tokens = line.strip().split()

                if(tokens[0] == 'L'):
                    tokens[1], tokens[3] = tokens[3], tokens[1]
                    if tokens[2] == "+" and tokens[4] == "+":
                        tokens[2] = "-"; tokens[4] = "-"
                    elif tokens[2] == "-" and tokens[4] == "-":
                        tokens[2] = "+"; tokens[4] = "+"
                    
                    print('\t'.join(tokens), file = ofile)
                
                print(line, file = ofile)

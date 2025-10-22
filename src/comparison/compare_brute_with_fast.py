import os
import json

data = ['MHC', 'C4-90', 'Mtb152m-p0a1', 'Mtb152m-p1a2', 'Mtb152p', 'human100', 'human100p10', 'Ecoli']

for i in range(len(data)):
    fast_path = f'/home/daanish/projects/panbubble/test/data/results/'
    out_path = f'/home/daanish/projects/Billi_data/results/small/Billi/{data[i]}/compare_with_vg.txt'
    
    id, rid = get_ids(gfa_file[i])

    lb = []
    with open(billi_path, 'r') as file:
        for line in file:
            tokens = line.strip().split()
            if tokens[0] == 'Total':
                continue
            g1 = tokens[0]
            g2 = tokens[2]
            id1 = id[g1]; id2 = id[g2]
            if id1 > id2:
                id1, id2 = id2, id1
                
            lb.append((id1, id2))     
    lb = sorted(lb)

    lv = []
    with open(vg_path, 'r') as file:
        for line in file:
            tokens = json.loads(line)
            id1 = int(tokens["start"]["node_id"]); id2 = int(tokens["end"]["node_id"])
            if id1 > id2:
                id1, id2 = id2, id1
                
            lv.append((id1, id2))

    lv = sorted(lv)

    # checking for duplicates
    for i in range(1, len(lb)):
        assert lb[i] != lb[i - 1]
    
    for i in range(1, len(lv)):
        assert lv[i] != lv[i - 1]

    with open(out_path, 'w') as ofile:
        print(f"vg count : {len(lv)}", file = ofile)
        print("Bubbles present in vg but not in Billi:", file = ofile)

        for i in range(len(lv)):
            if lv[i] not in lb:
                print(rid[lv[i][0]], rid[lv[i][1]], sep = ' ', file = ofile)

        print('====================================================', file = ofile)

        print(f"Billi count : {len(lb)}", file = ofile)
        print("Bubbles present in Billi but not in vg:", file = ofile)

        for i in range(len(lb)):
            if lb[i] not in lv:
                print(rid[lb[i][0]], rid[lb[i][1]], sep = ' ', file = ofile)
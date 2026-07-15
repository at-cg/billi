# Commands and Options

### decompose

Report panbubbles and hairpins in the input pangenome graph.

```
./billi decompose [OPTIONS] -i EC7.gfa > out.txt
```

*  `-i, --input FILE`:<br/>
   Input file in GFA format
*  `-e, --exact`:<br/>
   (optional) Use exact algorithm instead of the fast heuristic approach (default: heuristic)
*  `-m --minAlleles INT`:<br/>
    Minimum number of alleles required for reporting a panbubble or hairpin (default: 0)
*  `-w --iWalk`:<br/>
    Ignore haplotype walks (default: false)

### compact

Compute a compacted pangenome graph by collapsing non-branching linear paths.

```
./billi compact [OPTIONS] -i inputgraph.gfa -o compactgraph.gfa 
```

*  `-i, --input FILE`:<br/>
   Input file in GFA format
*  `-o, --output FILE`:<br/>
   Output file in GFA format
*  `-r --reverse`:<br/>
   Print edges in a reverse complement fashion as well (default:false)
*  `-n --numeric`:<br/>
   Modify node labels to numeric format (default: false)
*  `-s --self_loop`:<br/>
   Retain self loops (default: false)

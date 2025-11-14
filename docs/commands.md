# Commands and Options

1. [decompose](#decompose)

* `decompose`:<br/>
         For finding the panbubble boundaries given an input GFA file.

### decompose

```
billi decompose [OPTIONS] -i EC7.gfa -o EC7
```

*  `-i, --input FILE`:<br/>
   The input GFA file.
*  `-o, --output FILE`:<br/>
   The output directory where the results will be saved.
*  `-e, --exact [heuristic]`:<br/>
   (optional) Use the exact implementation to compute panbubbles.
*  `-c, --cycle-equivalent FLAG`:<br/>
   (optional) Whether cycle equivalent classes are to be reported?
*  `-r, --report-hairpins FLAG`:<br/>
   (optional) Whether hairpins are to be reported?

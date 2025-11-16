# Billi

<!-- ## **Dependencies**
> [pangene](https://github.com/lh3/pangene/tree/main) (For generating gfa file from the paf files)

> [CLI11.hpp](https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp) (copy in the ~/include directory) -->

## **Installation**

The binary will be saved in the bin folder.

```
git clone https://github.com/at-cg/billi.git
cd billi/src
make 
# This will download CLI11.hpp file in your ~/include directory if missing
```

## **Test run:**

### Enumerating panbubbles and hairpins
```
billi decompose -i <PATH_TO_GFA_FILE> -o <PATH_TO_THE_DIRECTORY_THAT_WILL_SAVE_THE_RESULTS> 
```
Check [here](docs/commands.md/#decompose) for the options. The output format is similar to [pangene](https://github.com/lh3/pangene/tree/main) but we only report the entrance edges along with the direction of traversal. 

There will be two files saved in the output directory:
```
summary.txt - This will list the summary of the input graph
panbubble.txt - This will list the panbubbles
```

### Example
```
billi decompose -i ../test/gfa_files/t2-1.gfa -c -r -o ../test/results/t2-1
```
<p align="center">
  <img src="docs/figures/t2-1.png" width="200">
  <br>
  <em>Figure 1: Bandage visualisation of t2-1 graph</em>
</p>

#### Billi output
```
<s6 <s4
>s1 >s3
```
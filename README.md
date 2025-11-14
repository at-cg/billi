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
Check [here](docs/commands.md/#decompose) for the options.

### Example
```
billi decompose -i ../test/gfa_files/t2-1.gfa -c -r -o ../test/results/t2-1
```
<p align="center">
  <img src="docs/figures/t2-1.png" width="400">
  <br>
  <em>Figure 1: Bandage visualisation</em>
</p>

#### Billi output
```
<s6 <s4
>s1 >s3
```
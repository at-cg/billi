# Panbubble

<!-- ## **Dependencies**
> [pangene](https://github.com/lh3/pangene/tree/main) (For generating gfa file from the paf files)

> [CLI11.hpp](https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp) (copy in the ~/include directory) -->

## **Installation**

The binary will be saved in the bin folder.

```
git clone https://github.com/at-cg/billi.git
cd billi/src
make (This will download CLI11.hpp file in your ~/include directory if missing)
```

## **Test run:**

### Enumerating panbubbles
```
billi decompose -i <PATH_TO_GFA_FILE> -o <PATH_TO_THE_DIRECTORY_THAT_WILL_SAVE_THE_RESULTS> 
```
Check [here](docs/commands.md/#decompose) for the options.

### Example
```
billi decompose -i ../test/gfa_files/EC7.gfa -c -r -o ../test/results/EC7
```
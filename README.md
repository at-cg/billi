# Panbubble

## **Dependencies**
> [pangene](https://github.com/lh3/pangene/tree/main) (For generating gfa file from the paf files)

> [CLI11.hpp](https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp) (copy in the ~/include directory)

## **Installation**
```
git clone https://github.com/at-cg/billi.git
cd billi/src
make
```

## **Test run:**

### Enumerating panbubbles
```
Assuming that you are in the src directory

./main decompose -i <PATH_TO_GFA_FILE> -o <PATH_TO_THE_DIRECTORY_THAT_WILL_SAVE_THE_RESULTS> 
```

## **Example:**
```
./main decompose -i ../test/gfa_files/EC7.gfa -c -r -o ../test/results/EC7
```
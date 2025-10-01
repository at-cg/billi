# Panbubble

## **Dependencies**
> [pangene](https://github.com/lh3/pangene/tree/main) (For generating gfa file from the paf files)

## **Installation**
```
git clone https://github.com/at-cg/Billi.git
cd src
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
./main decompose -i ../test/data/gfa_files/EC7.gfa -d 1 -c true -r true -p true -o ../test/data/results/EC7
```
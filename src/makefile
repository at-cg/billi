CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 -march=native -flto -funroll-loops -DNDEBUG -I ~/include# Compiler flags

TARGETS := main1 # Change this to your desired executable name 
SOURCES := main.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

main1: $(SOURCES)
		$(CXX) $(CXXFLAGS) -o $@ $^

clean:
		rm -f $(TARGETS)
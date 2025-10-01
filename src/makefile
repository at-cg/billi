CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 -I ~/include# Compiler flags

TARGETS := main # Change this to your desired executable name 
SOURCES := main.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

main1: $(SOURCES)
		$(CXX) $(CXXFLAGS) -o $@ $^

clean:
		rm -f $(TARGETS)
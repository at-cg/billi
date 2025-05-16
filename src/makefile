CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 # Compiler flags

TARGETS := main_opt # Change this to your desired executable name 
SOURCES := main_opt.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

main1: $(SOURCES)
		$(CXX) $(CXXFLAGS) -o $@ $^

clean:
		rm -f $(TARGETS)
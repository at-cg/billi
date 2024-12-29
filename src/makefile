CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 # Compiler flags

TARGETS := main # Change this to your desired executable name 
SOURCES := main.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

main: $(SOURCES)
		$(CXX) $(CXXFLAGS) -o $@ $^

clean:
		rm -f $(TARGETS)
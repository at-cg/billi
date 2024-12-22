CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 # Compiler flags
LDFLAGS := -L ~/lib # Library paths
INCLUDES := -I ~/include # Include paths

TARGETS := main # Change this to your desired executable name 
SOURCES := main.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

main: main.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		rm -f $(TARGETS)
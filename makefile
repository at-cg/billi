CXX := g++ # Compiler
# CXXFLAGS := -std=c++17 -g -O0 -O3 -I ~/include # GDB Compiler flags
CXXFLAGS := -std=c++17 -O3 -I ~/include # Compiler flags

TARGET := billi

SOURCES := src/main.cpp src/subcommand/decompose.cpp src/subcommand/compact.cpp src/include/commons.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
		$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp 
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -f $(TARGET) $(OBJECTS)
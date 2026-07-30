CXX := g++ # Compiler
DEBUG ?= 0
ifeq ($(DEBUG),1)
    # Use {make DEBUG=1} to enable this
    CXXFLAGS := -std=c++17 -g -O0 -I src/include
else
    # default compiler flags
    CXXFLAGS := -std=c++17 -O3 -I src/include
endif

LDFLAGS := -lstdc++fs

TARGET := billi

SOURCES := src/main.cpp src/subcommand/decompose.cpp src/subcommand/compact.cpp src/include/commons.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
		$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp 
		$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
		rm -f $(TARGET) $(OBJECTS)
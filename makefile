CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 -I ~/include # Compiler flags

TARGET_DIR := ..
TARGET := $(TARGET_DIR)/billi 
INCLUDE_DIR := ~/include

SOURCES := src/main.cpp

# CLI11_URL := https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp
# CLI11_HEADER := $(INCLUDE_DIR)/CLI11.hpp

# all: install_cli11 $(TARGET)

# install_cli11:
# 	@mkdir -p $(INCLUDE_DIR)
# 	@if [ ! -f $(CLI11_HEADER) ]; then \
# 	  echo "Downloading CLI11.hpp …"; \
# 	  curl -L $(CLI11_URL) -o $(CLI11_HEADER); \
# 	  echo "Saved to $(CLI11_HEADER)"; \
# 	else \
# 	  echo "CLI11.hpp already present"; \
# 	fi

$(TARGET): $(SOURCES)
# 		@mkdir -p $(TARGET_DIR)
		$(CXX) $(CXXFLAGS) -o $@ $^

clean:
		rm -f $(TARGET)

# Makefile (builds ./raytracer next to sources as to use as indicated in HW file)

CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall
LDFLAGS  :=

SRC := raytracer.cpp parser.cpp tinyxml2.cpp ppm.cpp
OBJ := $(SRC:.cpp=.o)
BIN := raytracer

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(BIN)

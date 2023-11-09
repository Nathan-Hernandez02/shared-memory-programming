CXX = g++
CXXFLAGS = -std=c++11 -Wall -fopenmp
SRC = numerical-integration.cpp
OUT = integration

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(OUT) $(SRC)

.PHONY: clean

clean:
	rm -f $(OUT)

run:
	./integration
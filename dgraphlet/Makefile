# Set to 1 to enable OpenMP
OPENMP_LIB ?= 1

all: d2s3 d2s3_b1 d2s3_b2 d3s3 d3s3_b1 d3s3_b2 d2s4

CXXFLAGS = -O3 -g -std=c++17
ifeq ($(OPENMP_LIB),1)
    CXXFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif
CXXFLAGS += -ffast-math -funroll-loops -fno-strict-aliasing \
            -fomit-frame-pointer -march=native
CXX = g++

# Target executables
d2s3: d2s3_b2.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d2s3 d2s3.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s3.cpp

d2s3_b1: d2s3_b1.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d2s3_b1 d2s3_b1.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s3.cpp

d2s3_b2: d2s3_b2.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d2s3_b2 d2s3_b2.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s3.cpp

d3s3: d3s3.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d3s3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d3s3 d3s3.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d3s3.cpp

d3s3_b1: d3s3_b1.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d3s3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d3s3_b1 d3s3_b1.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d3s3.cpp

d3s3_b2: d3s3_b2.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d3s3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d3s3_b2 d3s3_b2.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d3s3.cpp

d2s4: d2s4.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s4.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o d2s4 d2s4.cpp ./util/graph.cpp ./util/digraph.cpp ./util/counter_d2s4.cpp


clean:
	@echo "Removing executables"
	rm -f d2s3 d2s3_b1 d2s3_b2 d3s3 d3s3_b1 d3s3_b2 d2s4

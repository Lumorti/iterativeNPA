CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native
#CXXFLAGS=-g -fmax-errors=3 -Og -march=native -Wall
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/header_only_version/ 
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSOPTIM)
FILES=$(wildcard src/*.cpp)
ASOBJ=$(FILES:.cpp=.o)

all: run

run: $(ASOBJ)
	$(CXX) $(CXXFLAGS) -o run $(ASOBJ) $(LIBS)

src/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp -o src/main.o $(LIBSEIGEN) $(LIBSOPTIM)

src/mosek.o: src/mosek.cpp src/mosek.h
	$(CXX) $(CXXFLAGS) -c src/mosek.cpp -o src/mosek.o $(LIBSMOSEK) $(LIBSEIGEN)

src/optim.o: src/optim.cpp src/optim.h
	$(CXX) $(CXXFLAGS) -c src/optim.cpp -o src/optim.o $(LIBSOPTIM) $(LIBSEIGEN)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN)

clean:
	rm -f run src/*.o


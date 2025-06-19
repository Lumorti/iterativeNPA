CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native -fopenmp -pipe
DEBUGFLAGS=-g -fmax-errors=3 -Og -fopenmp
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/header_only_version/ 
LIBSSCS= -I${SCSHOME}/include/scs/ -L${SCSHOME}/lib/ -lscsdir
LIBSGUROBI= -I${GUROBIHOME}/include -L${GUROBIHOME}/lib -lgurobi_c++ -lgurobi110
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSOPTIM) $(LIBSSCS) $(LIBSGUROBI)
FILES=$(wildcard src/*.cpp)
ASOBJ=$(FILES:.cpp=.o)

all: run 

run: src/main.cpp ../PolyNC/builds/polyncNPA.o
	$(CXX) $(CXXFLAGS) -o run src/main.cpp ../PolyNC/builds/polyncNPA.o $(LIBS)

debug: src/main.cpp ../PolyNC/builds/polyncNPA.o
	$(CXX) $(DEBUGFLAGS) -o run src/main.cpp ../PolyNC/builds/polyncNPA.o $(LIBS)

clean:
	rm -f run src/*.o


CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native -fopenmp
#CXXFLAGS=-g -fmax-errors=3 -O0 -fopenmp
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/header_only_version/ 
LIBSSPECTRA= -I${SPECTRAHOME}/include/
LIBSSCS= -I${SCSHOME}/include/scs/ -L${SCSHOME}/lib/ -lscsdir
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSOPTIM) $(LIBSSCS) $(LIBSSPECTRA)
FILES=$(wildcard src/*.cpp)
ASOBJ=$(FILES:.cpp=.o)

all: run

run: $(ASOBJ)
	$(CXX) $(CXXFLAGS) -o run $(ASOBJ) $(LIBS)

src/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp -o src/main.o $(LIBSEIGEN) $(LIBSOPTIM) $(LIBSSPECTRA)

src/optMOSEK.o: src/optMOSEK.cpp src/optMOSEK.h
	$(CXX) $(CXXFLAGS) -c src/optMOSEK.cpp -o src/optMOSEK.o $(LIBSMOSEK) $(LIBSEIGEN)

src/optOptim.o: src/optOptim.cpp src/optOptim.h
	$(CXX) $(CXXFLAGS) -c src/optOptim.cpp -o src/optOptim.o $(LIBSOPTIM) $(LIBSEIGEN)

src/optSCS.o: src/optSCS.cpp src/optSCS.h
	$(CXX) $(CXXFLAGS) -c src/optSCS.cpp -o src/optSCS.o $(LIBSSCS) $(LIBSEIGEN)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN)

clean:
	rm -f run src/*.o


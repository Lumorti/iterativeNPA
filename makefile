CXX=g++
CXXFLAGS=-Wfatal-errors -O3 -march=native
DEBUGFLAGS=-Wall -g -pg -Wfatal-errors
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/include/ -L${OPTIMHOME} -loptim
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSOPTIM)

all: 
	$(CXX) $(CXXFLAGS) -o run main.cpp $(LIBS)

debug:
	$(CXX) $(DEBUGFLAGS) -o run main.cpp $(LIBS)


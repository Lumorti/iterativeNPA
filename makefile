CXX=g++
CXXFLAGS=-Wfatal-errors -O3 -march=native
DEBUGFLAGS=-Wall -g -pg -Wfatal-errors
LIBSNN= -I/home/luke/Documentos/MiniDNN/include
LIBSEIGEN= -I/usr/include/eigen3
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSNN)

all: 
	$(CXX) $(CXXFLAGS) -o run main.cpp $(LIBS)

debug:
	$(CXX) $(DEBUGFLAGS) -o run main.cpp $(LIBS)


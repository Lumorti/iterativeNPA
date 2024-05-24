CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native
DEBUGFLAGS=-g -fmax-errors=3 -O0
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/header_only_version/ 
LIBSAUTODIFF= -I${AUTODIFFHOME}/
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSOPTIM) $(LIBSAUTODIFF)

all: 
	$(CXX) $(CXXFLAGS) -o run main.cpp $(LIBS)

debug:
	$(CXX) $(DEBUGFLAGS) -o run main.cpp $(LIBS)


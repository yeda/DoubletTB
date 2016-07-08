CONFIG=root-config
CXXFLAGS=$(shell $(CONFIG) --cflags)
LIBS=$(shell $(CONFIG) --glibs)
LDFLAGS=$(shell $(CONFIG) --ldflags)
CXX=g++
# enable debug settings for GCC
ADDCXXFLAGS=-ggdb -O0 -std=c++0x
#ADDCXXFLAGS=-ggdb -O0 -std=c++11
# OR enable performance optimizations
#ADDCXXFLAGS=-g -O2
VPATH=./
HDIR=./

HDRS= $(HDIR)/Settings.h $(HDIR)/HitMaker.h $(HDIR)/Alignment.h $(HDIR)/Analyze.h 

RUNALLHDRS= $(HDIR)/RunAll.h

HITMAKEROBJS=	HitMaker.o VecVecDict.o 
ALIGNMENTOBJS= Alignment.o VecVecDict.o
ANALYZEOBJ= Analyze.o VecVecDict.o
RUNALLOBJ= RunAll.o

all: hitmaker alignment analyze runall

runall: $(RUNALLOBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(ADDCXXFLAGS) $(RUNALLOBJ) $(LDFLAGS) $(LIBS)

hitmaker: $(HITMAKEROBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(ADDCXXFLAGS) $(HITMAKEROBJS) $(LDFLAGS) $(LIBS)

alignment: $(ALIGNMENTOBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(ADDCXXFLAGS) $(ALIGNMENTOBJS) $(LDFLAGS) $(LIBS)

analyze: $(ANALYZEOBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(ADDCXXFLAGS) $(ANALYZEOBJ) $(LDFLAGS) $(LIBS)

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) -c $< 

clean:
	rm -rf hitmaker alignment analyze runall $(RUNALLOBJ) $(ANALYZEOBJ) $(ALIGNMENTOBJS) $(HITMAKEROBJS) $(OBJS) 

SRCS = $(wildcard *.cc)
TRGS = $(patsubst %.cc, ./%, $(SRCS))

ROOTLIBS=`root-config --cflags --glibs`
CXXFLAGS=--std=c++17

all: $(TRGS)

./% : %.cc
	g++ -g -o $@ $(CXXFLAGS) $(LIBS) $(ROOTLIBS) $(INCDIR) $<


include ../Makefile.include

LFLAGS  += -DEXPANSION=40
LFLAGS 	+= -DIneJ
LFLAGS  += -DFP64

serial	: main.cxx
	$(CXX) $? $(LFLAGS)
	./a.out -n 100000

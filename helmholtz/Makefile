include ../Makefile.include

LFLAGS	+= -DFP64

serial	: main.cxx
	$(CXX) $? $(LFLAGS) -DEXPANSION=40
	./a.out -n 100000

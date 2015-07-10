include ../Makefile.include

LFLAGS	+= -DFP64 -DSpherical -DEXPANSION=10

serial	: main.cxx
	$(CXX) $? $(LFLAGS)
	./a.out -n 100000

include ../Makefile.include

LFLAGS	+= -DFP64

serial	: main.cxx
	$(CXX) $? $(LFLAGS)
	./a.out -n 100000

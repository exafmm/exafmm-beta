include ../Makefile.include

LFLAGS 	+= -DIneJ

serial	: main.cxx
	$(CXX) $? $(LFLAGS)
	./a.out -n 100000

include ../Makefile.include
LFLAGS	+= -DSTACK
#LFLAGS	+= -DIneJ
#LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

a	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	./a.out

b	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

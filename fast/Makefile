include ../Makefile.include
#LFLAGS	+= -DIneJ
#LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

a	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DDIRECT
	./a.out

b	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT
	./a.out

d	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

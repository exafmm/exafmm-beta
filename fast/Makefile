include ../Makefile.include
LFLAGS	+= -DCOMkernel
LFLAGS	+= -DSTACK
#LFLAGS	+= -DIneJ
LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

a	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	MYTH_WORKER_NUM=1 ./a.out

b	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY
	MYTH_WORKER_NUM=1 ./a.out

c	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

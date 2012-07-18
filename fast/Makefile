include ../Makefile.include
LFLAGS	+= -DCOMkernel
LFLAGS	+= -DERROR_OPT
LFLAGS	+= -DUSE_RMAX
LFLAGS	+= -DUSE_BMAX
LFLAGS	+= -DCOMcenter
LFLAGS	+= -DDUAL
LFLAGS	+= -DIneJ
#LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

a	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	MYTH_WORKER_NUM=12 mpirun -np 2 ./a.out

b	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY
	MYTH_WORKER_NUM=12 ./a.out

c	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

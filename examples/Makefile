include ../Makefile.include
LFLAGS	+= -DCOMkernel
LFLAGS	+= -DERROR_OPT
LFLAGS	+= -DUSE_RMAX
LFLAGS	+= -DUSE_BMAX
LFLAGS	+= -DCOMcenter
LFLAGS	+= -DDUAL
LFLAGS	+= -DIneJ
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

serial	: serial.cxx
	make clean
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	MYTH_WORKER_NUM=12 ./a.out

parallel: parallel.cxx
	make clean
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	MYTH_WORKER_NUM=12 mpirun -np 4 ./a.out

b	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY
	MYTH_WORKER_NUM=12 mpirun -np 4 ./a.out

include ../Makefile.include
#LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

a	: main.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) $(VFLAGS)
#	mpirun -np 2 -machinefile ~/hosts ./a.out
	mpirun -np 1 ./a.out

b	: main.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

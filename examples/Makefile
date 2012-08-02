include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	MYTH_WORKER_NUM=12 ./a.out

parallel: parallel.cxx
	make clean
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	MYTH_WORKER_NUM=12 mpirun -np 16 ./a.out

wrapper: wrapper.cxx
	make clean
	make -C ../wrapper libcoulomb.a
	$(CXX) $? $(LFLAGS) -lcoulomb
	MYTH_WORKER_NUM=12 mpirun -np 4 ./a.out
	$(RM) libcoulomb.a
	make -C ../wrapper libcoulomb.so
	$(CXX) $? $(LFLAGS) -lcoulomb
	MYTH_WORKER_NUM=12 mpirun -np 4 ./a.out

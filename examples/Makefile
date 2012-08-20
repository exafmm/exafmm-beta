include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	./a.out

serial_argv : serial_argv.cxx
	make clean
	$(CXX) -DCMDLINE_ARGS $? $(LFLAGS) $(VFLAGS)
	./a.out --mutual 1 --parallelEverything 1 --splitBothThreshold 1000 --ncrit 10

parallel: parallel.cxx
	make clean
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 32 ./a.out

wrapper: wrapper.cxx
	make clean
	make -C ../wrapper libcoulomb.a
	$(CXX) $? $(LFLAGS) -lcoulomb
	mpirun -np 4 ./a.out

include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	./a.out

serial_argv : serial_argv.cxx
	make clean
	$(CXX) -DCMDLINE_ARGS $? $(LFLAGS)
	./a.out --mutual 0 --parallelEverything 1 --splitBothThreshold 1000 --ncrit 10 --simdize sse 

parallel: parallel.cxx
	make clean
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 32 ./a.out

wrapper: wrapper.cxx
	make clean
	make -C ../wrapper libcoulomb.a
	$(CXX) $? $(LFLAGS) -lcoulomb
	mpirun -np 4 ./a.out

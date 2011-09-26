include ../Makefile

main: main.cxx $(OBJECT)
	rm statistics.dat
	@$(CXX) $? $(LFLAGS) -lfftw3
	mpirun -np 256 -machinefile ~/hosts ./a.out

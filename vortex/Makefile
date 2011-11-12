include ../Makefile

main: main.cxx $(OBJECT)
	rm -f statistics.dat
	@$(CXX) $? $(LFLAGS) -lfftw3
	mpirun -np 4 ./a.out
#	mpirun -np 16 -machinefile ~/hosts ./a.out

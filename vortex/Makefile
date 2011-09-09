include ../Makefile

main: main.cxx $(OBJECT)
	rm -f statistics.dat
	@$(CXX) $? $(LFLAGS) -lfftw3
	mpirun -np 1 ./a.out

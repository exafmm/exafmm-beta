include ../Makefile

main: main.cxx $(OBJECT)
	rm -f statistics.dat
	@$(CXX) $? $(LFLAGS) -lfftw3
	./a.out

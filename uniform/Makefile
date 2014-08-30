include ../Makefile.include

serial	: main.cxx
	$(CXX) $? $(LFLAGS) -DSerial
	./a.out -n 100000 -c 100 -i 0

parallel: main.cxx
	$(CXX) $? $(LFLAGS)
	mpirun -np 2 ./a.out -n 100000 -c 100 -i 0

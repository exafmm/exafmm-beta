include ../Makefile

coulombVdW: coulombVdW.o $(OBJECT)
	ar ruv libcoulombVdW.a $?
	ranlib libcoulombVdW.a

testlib: testlib.cxx
	make cleanall
	rm -f libcoulombVdW.a
	make coulombVdW
	@$(CXX) $? libcoulombVdW.a $(LFLAGS)
	mpirun -np 2 ./a.out

ewald: ewald.cxx mr3.o mr3_host.o
	$(CXX) $? $(LFLAGS)
	./a.out 0 12
	./a.out 1 12

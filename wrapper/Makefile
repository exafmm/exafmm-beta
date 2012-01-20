include ../Makefile.include

ifeq ($(DEVICE),gpu)
MR3 = mr3.o
else
MR3 = mr3_host.o
endif

coulomb: coulomb.o $(OBJECT)
	ar ruv libcoulomb.a $?
	ranlib libcoulomb.a

coulombVdW: coulombVdW.o $(OBJECT)
	ar ruv libcoulombVdW.a $?
	ranlib libcoulombVdW.a

test_coulomb: test_coulomb.cxx $(MR3)
	rm -f libcoulomb.a
	make coulomb
	$(CXX) $? libcoulomb.a $(LFLAGS)
	mpirun -np 2 ./a.out

test_coulombVdW: test_coulombVdW.cxx $(MR3)
	rm -f libcoulombVdW.a
	make coulombVdW
	$(CXX) $? libcoulombVdW.a $(LFLAGS)
	mpirun -np 2 ./a.out

ewald: ewald.cxx $(MR3)
	rm -f libcoulombVdW.a
	make coulombVdW
	$(CXX) $? libcoulombVdW.a $(LFLAGS)
	./a.out 0 12
	./a.out 1 12


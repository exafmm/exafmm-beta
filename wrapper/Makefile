include ../Makefile.include

ifeq ($(DEVICE),gpu)
MR3 = mr3.o
else
MR3 = mr3_host.o
endif

libcoulomb.a: coulomb.o $(OBJECT)
	ar ruv libcoulomb.a $?
	ranlib libcoulomb.a

libcoulomb.so: coulomb.o $(OBJECT)
	$(CXX) $? -shared -o libcoulomb.so $(LFLAGS)

libcoulombVdW.a: coulombVdW.o $(OBJECT)
	ar ruv libcoulombVdW.a $?
	ranlib libcoulombVdW.a

libcoulombVdW.so: coulombVdW.o $(OBJECT)
	$(CXX) $? -shared -o libcoulombVdW.so $(LFLAGS)

test_coulomb: test_coulomb.cxx $(MR3)
	rm -f libcoulomb.*
	make libcoulomb.a
	$(CXX) $? -L. -lcoulomb $(LFLAGS)
	mpirun -np 2 ./a.out
	rm -f libcoulomb.*
	make libcoulomb.so
	$(CXX) $? -L. -lcoulomb $(LFLAGS)
	mpirun -np 2 ./a.out

test_coulombVdW: test_coulombVdW.cxx $(MR3)
	rm -f libcoulombVdW.*
	make libcoulombVdW.a
	$(CXX) $? -L. -lcoulombVdW $(LFLAGS)
	mpirun -np 2 ./a.out
	rm -f libcoulombVdW.*
	make libcoulombVdW.so
	$(CXX) $? -L. -lcoulombVdW $(LFLAGS)
	mpirun -np 2 ./a.out

ewald: ewald.cxx $(MR3)
	rm -f libcoulombVdW.a
	make coulombVdW
	$(CXX) $? libcoulombVdW.a $(LFLAGS)
	./a.out 0 12
	./a.out 1 12

include ../Makefile

ifeq ($(DEVICE),gpu)
MR3 = mr3.o
else
MR3 = mr3_host.o
endif

coulombVdW: coulombVdW.o $(OBJECT)
	ar ruv libcoulombVdW.a $?
	ranlib libcoulombVdW.a

testlib: testlib.cxx $(MR3)
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


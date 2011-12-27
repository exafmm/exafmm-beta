include ../Makefile

coulombVdW: coulombVdW.o $(OBJECT)
	ar ruv libcoulombVdW.a $?
	ranlib libcoulombVdW.a

testlib: testlib.cxx
	make cleanall
	rm -f libcoulombVdW.a
	make coulombVdW
	@$(CXX) $? libcoulombVdW.a $(LFLAGS)
	./a.out

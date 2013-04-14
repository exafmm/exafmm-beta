include ../Makefile.include

libcoulomb.a: coulomb.o
	ar ruv libcoulomb.a $?
	ranlib libcoulomb.a

libcoulomb.so: coulomb.o
	$(CXX) $? -shared -o libcoulomb.so $(LFLAGS)

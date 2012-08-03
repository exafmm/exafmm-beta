include ../Makefile.include

libcoulomb.a: coulomb.o
	ar ruv libcoulomb.a $?
	ranlib libcoulomb.a
	mkdir -p $(EXAFMM_LIBRARY_PATH)
	mv libcoulomb.a $(EXAFMM_LIBRARY_PATH)

libcoulomb.so: coulomb.o
	$(CXX) $? -shared -o libcoulomb.so $(LFLAGS)
	mkdir -p $(EXAFMM_LIBRARY_PATH)
	mv libcoulomb.so $(EXAFMM_LIBRARY_PATH)

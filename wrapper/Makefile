include ../Makefile.include

ifeq ($(DEVICE),GPU)
MR3 = mr3.o
else
MR3 = mr3_host.o
endif

libcoulomb.a: coulomb.o
	ar ruv libcoulomb.a $?
	ranlib libcoulomb.a
	mv libcoulomb.a $(EXAFMM_LIBRARY_PATH)

libcoulomb.so: coulomb.o
	$(CXX) $? -shared -o libcoulomb.so $(LFLAGS)
	mv libcoulomb.so $(EXAFMM_LIBRARY_PATH)

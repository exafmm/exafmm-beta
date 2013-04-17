include ../Makefile.include

libfmm.a: fmm.o
	ar ruv libfmm.a $?
	ranlib libfmm.a

libfmm.so: fmm.o
	$(CXX) $? -shared -o libfmm.so $(LFLAGS)

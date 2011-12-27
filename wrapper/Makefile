include ../Makefile

host: host.o ../kernel/cpu$(EXPAND)Laplace.o ../kernel/cpuCoulombVdW.o
	rm -f libFMM.a
	ar ruv libFMM.a $?
	ranlib libFMM.a

device: device.o ../kernel/gpu$(EXPAND)Laplace.o ../kernel/gpuCoulombVdW.o
	rm -f libFMM.a
	ar ruv libFMM.a $?
	ranlib libFMM.a

testlib: testlib.cxx
	make cleanall
	make host
	@$(CXX) $? libFMM.a $(LFLAGS)
	./a.out

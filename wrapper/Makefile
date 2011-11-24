include ../Makefile

host: host.o ../kernel/cpu$(EXPAND)Laplace.o ../kernel/cpu$(EXPAND)CoulombVdW.o
	ar ruv libFMM.a $?
	ranlib libFMM.a

device: device.o ../kernel/gpu$(EXPAND)Laplace.o ../kernel/gpu$(EXPAND)CoulombVdW.o
	ar ruv libFMM.a $?
	ranlib libFMM.a

testlib: testlib.cxx
	make host
	@$(CXX) $? libFMM.a
	./a.out

include ../Makefile

static: host.o $(OBJECT)
	ar ruv libFMM.a $?
	ranlib libFMM.a

nbody: nbody.cxx
	make static
	@$(CXX) $? libFMM.a
	./a.out

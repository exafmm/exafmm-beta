include ../Makefile

host: host.o $(OBJECT)
	ar ruv libFMM.a $?
	ranlib libFMM.a

nbody: nbody.cxx
	make host
	@$(CXX) $? libFMM.a
	./a.out

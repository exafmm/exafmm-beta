include ../Makefile

static: Laplace.o $(OBJECT)
	ar ruv ../lib/lib$(DEVICE)Laplace.a $?
	ranlib ../lib/lib$(DEVICE)Laplace.a

nbody: nbody.cxx
	@$(CXX) $? ../lib/libcpuLaplace.a
	./a.out

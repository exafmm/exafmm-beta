include ../Makefile

nbody: nbody.cxx
	@$(CXX) $? ../lib/libcpuLaplace.a $(LFLAGS)
	./a.out

vortex: vortex.cxx
	@$(CXX) $? ../lib/libcpuBiotSavart.a ../lib/libcpuStretching.a $(LFLAGS)
	./a.out

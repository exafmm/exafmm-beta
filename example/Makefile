include ../Makefile

nbody: nbody.cxx
	@$(CXX) $? ../lib/libgpuLaplace.a $(LFLAGS)
	./a.out

vortex: vortex.cxx
	@$(CXX) $? ../lib/libcpuBiotSavart.a ../lib/libcpuStretching.a ../lib/libcpuGaussian.a $(LFLAGS)
	./a.out

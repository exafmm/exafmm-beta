include ../Makefile

nbody: nbody.cxx
	@$(CXX) $? ../lib/libgpuLaplace.a $(LFLAGS)
	./a.out

vortex: vortex.cxx
	@$(CXX) $? ../lib/libgpuBiotSavart.a ../lib/libgpuStretching.a ../lib/libgpuGaussian.a $(LFLAGS)
	./a.out

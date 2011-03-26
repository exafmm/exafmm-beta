include ../Makefile

nbody: nbody.cxx
	@$(CXX) $? ../lib/libcpuLaplace.a
	./a.out

vortex: vortex.cxx
	@$(CXX) $? ../lib/libgpuBiotSavart.a ../lib/libgpuStretching.a -lcudart
	./a.out

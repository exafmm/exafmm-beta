include ../Makefile

nbody: nbody.cxx
	@$(CXX) $? ../lib/libcpuLaplace.a
	./a.out

vortex: vortex.cxx
	@$(CXX) $? ../lib/libcpuBiotSavart.a ../lib/libcpuStretching.a -lcudart
	./a.out

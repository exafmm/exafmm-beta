include ../Makefile

nbody: nbody.cxx
	@$(CXX) $? ../lib/libcpuLaplace.a
	./a.out

vortex: vortex.cxx
	@$(CXX) $? -L../lib -lcpuStretching -lcpuBiotSavart
	./a.out

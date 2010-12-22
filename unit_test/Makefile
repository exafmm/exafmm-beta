include ../Makefile

link    : link.cxx
	@$(CXX) $<
	./a.out

vtk     : vtk.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

morton  : morton.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	valgrind ./a.out

topdown : topdown.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	valgrind ./a.out

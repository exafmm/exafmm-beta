include ../Makefile

link    : link.cxx
	@$(CXX) -g -pg $< $(OFLAGS) $(LFLAGS)
	./a.out
	gprof -b a.out gmon.out

vtk     : vtk.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	./a.out

bottomup: bottomup.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

topdown : topdown.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

mpi     : mpi.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 128 ./a.out

alltoall: alltoall.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 8 ./a.out

alltoallv: alltoallv.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 4 ./a.out

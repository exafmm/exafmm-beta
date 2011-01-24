include ../Makefile

morton  : morton.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	./a.out

bottomup : bottomup.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

topdown  : topdown.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

kernel   : kernel.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

approx   : approx.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

direct   : direct.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

taylor   : taylor.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

vtk      : vtk.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	./a.out

mpi      : mpi.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) -g -pg
	mpirun -np 128 ./a.out
	gprof -b a.out gmon.out

alltoall : alltoall.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 8 ./a.out

alltoallv: alltoallv.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 4 ./a.out

median   : median.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 8 ./a.out

all      :
	make morton
	make bottomup
	make topdown
	make kernel
	make approx
	make direct
	make taylor
	make vtk
	make mpi
	make alltoall
	make alltoallv
	make median

include ../Makefile

index: index.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	./a.out

bottomup: bottomup.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

topdown: topdown.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

kernel: kernel.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

approx: approx.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

direct: direct.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

taylor: taylor.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	./a.out

vtk: vtk.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	./a.out

mpi: mpi.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) -g -pg
	mpirun -np 128 ./a.out
	gprof -b a.out gmon.out

shift: shift.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	mpirun -np 4 ./a.out

alltoall: alltoall.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 8 ./a.out

alltoallv: alltoallv.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS)
	mpirun -np 4 ./a.out

nth_element: nth_element.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	mpirun -np 8 ./a.out

multisection: multisection.cxx
	@$(CXX) $< $(OFLAGS) $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

all      :
	make index
	make bottomup
	make topdown
	make kernel
	make approx
	make direct
	make taylor
	make vtk
	make mpi
	make shift
	make alltoall
	make alltoallv
	make nth_element
	make multisection

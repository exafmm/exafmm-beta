include ../Makefile

index: index.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	./a.out

topdown: construct.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) -DTOPDOWN
	./a.out

bottomup: construct.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS)
	./a.out

kernel: kernel.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS)
	./a.out

serialrun: serialrun.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) -DTOPDOWN
	./a.out

direct_gpu: direct_gpu.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS)
	./a.out

taylor: taylor.cxx
	@$(CXX) $? $(LFLAGS)
	./a.out

vtk: vtk.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	./a.out

mpi: mpi.cxx
	@$(CXX) $? $(LFLAGS) -g -pg
	mpirun -np 128 ./a.out
	gprof -b a.out gmon.out

shift: shift.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 4 ./a.out

nth_element: nth_element.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

bisection: bisection.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

let: let.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

parallelrun: parallelrun.cxx $(KERNEL)
	@$(CXX) $? $(LFLAGS)
	mpirun -np 7 ./a.out

all      :
	make index
	make topdown
	make bottomup
	make kernel
	make serialrun
	make direct_gpu
	make taylor
	make vtk
	make mpi
	make shift
	make nth_element
	make bisection
	make let
	make parallelrun

include ../Makefile

index: index.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	./a.out

topdown: construct.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) -DTOPDOWN
	./a.out

bottomup: construct.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS)
	./a.out

kernel: kernel.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS)
	./a.out

serialrun: serialrun.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) -DTOPDOWN
	./a.out

unsort: unsort.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) -DTOPDOWN
	./a.out

Nscaling: Nscaling.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS)
	./a.out

cuprintf: cuprintf.cu
	@$(NVCC) $? $(LFLAGS)
	./a.out

direct_gpu: direct_gpu.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS)
	./a.out

mpi: mpi.cxx
	@$(CXX) $? $(LFLAGS) -g -pg
	mpirun -np 128 ./a.out
	gprof -b a.out gmon.out

check_gpus: check_gpus.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS)
#	mpirun -np 4 ./a.out

shift: shift.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 4 ./a.out

nth_element: nth_element.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

bisection: bisection.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

let: let.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 7 ./a.out

parallelrun: parallelrun.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS)
#	mpirun -np 4 ./a.out

vtk: vtk.cxx $(OBJECT)
	@$(CXX) $? $(LFLAGS) $(VFLAGS)
	./a.out

wrapper: wrapper.o $(OBJECT)
	ar ruv ../lib/lib$(DEVICE)$(KERNEL).a $?
	ranlib ../lib/lib$(DEVICE)$(KERNEL).a

all      :
	make index
	make topdown
	make bottomup
	make kernel
	make serialrun
	make direct_gpu
	make mpi
	make shift
	make nth_element
	make bisection
	make let
	make parallelrun

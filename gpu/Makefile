CXX = g++ -fopenmp -O3
NVCC = nvcc -g -arch=sm_35 -use_fast_math -Iinclude -Xcompiler "-fopenmp -O3"

ORG = approximate.o \
	main.o \
	scanKernels.o \
	tree.o \
	thrust.o

.SUFFIXES: .o .cpp .cu

all: $(ORG)
	$(CXX) $^ -lcudart
	./a.out

%.o: %.cu
	$(NVCC) -c $< -o $@

%.th_o: %.cu
	$(NVCC) -c $< -o $@

$(ORG): include/*.h

clean:
	rm -f approximate.o main.o scanKernels.o tree.o
cleanall:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
commit:
	@make -C .. commit
save:
	@make -C .. save
revert:
	@make -C .. revert

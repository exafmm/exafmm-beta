CUDA_INSTALL_PATH  = /usr/local/cuda

CXX = g++ -fopenmp -O3
NVCC = nvcc --compiler-bindir=/usr/bin/g++-4.4 -arch sm_21 -use_fast_math -Iinclude -Xcompiler "-fopenmp -O3"

ORG = approximate.o \
	main.o \
	scanKernels.o \
	tree.o
TORG = sortKernels.th_o

.SUFFIXES: .o .cpp .cu

all: $(ORG) $(TORG)
	$(CXX) $^ -lcudart -L$(CUDA_INSTALL_PATH)/lib64
	./a.out

%.o: %.cu
	$(NVCC) -c $< -o $@

%.th_o: %.cu
	$(NVCC) -c $< -o $@

$(ORG): include/*.h

clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
cleanall:
	make clean
	rm *.th_o
commit:
	@make -C .. commit
save:
	@make -C .. save
revert:
	@make -C .. revert

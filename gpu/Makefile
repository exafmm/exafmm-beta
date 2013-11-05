CXX = g++ -fopenmp -O3
NVCC = nvcc -g -arch=sm_35 -use_fast_math -Xcompiler "-fopenmp -O3"

ORG = approximate.o \
	serial.o \
	scanKernels.o \
	tree.o \
	thrust.o

.SUFFIXES: .o .cpp .cu

serial: $(ORG)
	$(CXX) $^ -lcudart
	./a.out

%.o: %.cu
	$(NVCC) -c $< -o $@

%.th_o: %.cu
	$(NVCC) -c $< -o $@

clean:
	rm -f approximate.o serial.o scanKernels.o tree.o
cleanall:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
commit:
	@make -C .. commit
save:
	@make -C .. save
revert:
	@make -C .. revert

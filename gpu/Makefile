.SUFFIXES: .cu

NVCC = /usr/local/cuda-5.5/bin/nvcc -g -use_fast_math -arch=sm_35 -rdc=true
LFLAGS = -lcudadevrt
LFLAGS += -DMASS # Use all positive sources

SRC = thrust.cu serial.cu

OBJ = $(SRC:%.cu=%.o)

# Test for serial FMM
serial: $(OBJ)
	$(NVCC) $^ $(LFLAGS)
	./a.out

%.o: %.cu
	$(NVCC) -c $< $(LFLAGS)

clean:
	rm -f serial.o

cleanall:
	rm -f *.o *.out

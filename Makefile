.SUFFIXES: .cxx .cu .o

DEVICE  = cpu
#DEVICE  = gpu

EXPAND  = Cartesian
#EXPAND  = Spherical

KERNEL  = Laplace

CXX     = mpicxx -mpreferred-stack-boundary=4 -ggdb3 -Wall -Wextra -Winit-self -Wshadow -O2 -fPIC -fopenmp\
	-ffast-math -funroll-loops -fforce-addr -rdynamic -D_FILE_OFFSET_BITS=64\
	-I../include -I/usr/include/vtk-5.2 -L/usr/lib/vtk-5.2
NVCC    = nvcc --ptxas-options=-v -O3 -use_fast_math -arch=sm_11\
	-I../include -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc
LFLAGS  = -L$(CUDA_INSTALL_PATH)/lib64 -L$(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lcutil_x86_64 -lstdc++ -ldl -lm\
	-D$(DEVICE) -D$(EXPAND) -D$(KERNEL)
VFLAGS  = -lvtkHybridTCL -lvtkWidgetsTCL -DVTK
OBJECT  = ../kernel/$(DEVICE)$(EXPAND)$(KERNEL).o ../kernel/$(DEVICE)Evaluator.o

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)
.cu.o   :
	$(NVCC) -c $? -o $@ $(LFLAGS)
cleanall:
	rm -rf *.o *.out *.sum ../kernel/*.o
save    :
	make cleanall
	tar zcvf ../../exafmm.tgz ../../exafmm

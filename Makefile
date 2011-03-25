.SUFFIXES: .cxx .cu .o

#DEVICE  = cpu
DEVICE  = gpu

#EXPAND  = Cartesian
EXPAND  = Spherical

KERNEL  = Laplace
#KERNEL  = BiotSavart
#KERNEL  = Stretching

CXX     = mpicxx -mpreferred-stack-boundary=4 -ggdb3 -Wall -Wextra -Winit-self -Wshadow -O2 -fPIC -fopenmp\
	-ffast-math -funroll-loops -fforce-addr -rdynamic -D_FILE_OFFSET_BITS=64\
	-I../include -I/usr/include/vtk-5.2 -L/usr/lib/vtk-5.2
NVCC    = nvcc --ptxas-options=-v -O3 -use_fast_math -arch=sm_11\
	-I../include -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc
LFLAGS  = -L$(CUDA_INSTALL_PATH)/lib64 -L$(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lcutil_x86_64 -lstdc++ -ldl -lm\
	-D$(DEVICE) -D$(EXPAND) -D$(KERNEL)
VFLAGS  = -lvtkHybridTCL -lvtkWidgetsTCL -DVTK
OBJECT  = ../kernel/$(DEVICE)$(EXPAND)$(KERNEL).o ../kernel/$(DEVICE)Evaluator.o
ifeq ($(DEVICE),gpu)
OBJECT  += ../kernel/$(KERNEL)Evaluator.o
endif

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)
.cu.o   :
	$(NVCC) -c $? -o $@ $(LFLAGS)
cleanall:
	rm -f ../unit_test/*.o ../unit_test/*.out ../unit_test/*.sum ../unit_test/time
	rm -f ../example/*.o ../example/*.out ../example/time ../kernel/*.o ../wrapper/*.o
cleanlib:
	rm -f ../lib/*.a ../lib/*.so
save    :
	make cleanall
	make cleanlib
	tar zcvf ../../exafmm.tgz ../../exafmm

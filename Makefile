.SUFFIXES: .cxx .cu .o

DEVICE  = cpu
#DEVICE  = gpu

EXPAND  = Cartesian
#EXPAND  = Spherical

CXX     = mpicxx -mpreferred-stack-boundary=4 -ggdb3 -Wall -Wextra -Winit-self -Wshadow -O3 -fPIC -fopenmp\
	-ffast-math -funroll-loops -fforce-addr -rdynamic -D_FILE_OFFSET_BITS=64\
	-I. -I../include -I/usr/local/fftw/include -I/usr/include/vtk-5.2 -L/usr/lib/vtk-5.2
#CXX     = mpicxx -O2 -fPIC -openmp -I../include -I/usr/local/fftw/include -I/usr/include/vtk-5.2 -L/usr/lib/vtk-5.2
NVCC    = nvcc -Xcompiler -fopenmp --ptxas-options=-v -O3 -use_fast_math -arch=sm_13\
	-I../include -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc
LFLAGS  = -L$(CUDA_INSTALL_PATH)/lib64 -L$(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lcutil_x86_64 -lstdc++ -ldl -lm\
	-D$(DEVICE) -D$(EXPAND)
#VFLAGS  = -lvtkHybridTCL -lvtkWidgetsTCL -DVTK
OBJECT  = ../kernel/$(DEVICE)$(EXPAND)Laplace.o ../kernel/$(DEVICE)$(EXPAND)BiotSavart.o\
	../kernel/$(DEVICE)$(EXPAND)Stretching.o ../kernel/$(DEVICE)$(EXPAND)Gaussian.o\
	../kernel/$(DEVICE)$(EXPAND)CoulombVdW.o

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)
.cu.o   :
	$(NVCC) -c $? -o $@ $(LFLAGS)
cleanall:
	rm -f ../unit_test/*.o ../unit_test/*.out ../unit_test/*.sum ../unit_test/time ../unit_test/direct0*
	rm -f ../example/*.o ../example/*.out ../example/time ../kernel/*.o ../wrapper/*.o
	rm -f ../wrapper/*.o ../wrapper/*.out ../wrapper/*.a ../wrapper/time
	rm -f ../vortex/*.o ../vortex/*.out ../vortex/time
	rm -f ../fast/*.o ../fast/*.out ../fast/time
save    :
	make cleanall
	tar zcvf ../../exafmm.tgz ../../exafmm

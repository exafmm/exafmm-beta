.SUFFIXES: .cxx .cu .o
.PHONY: docs

CUDA_INSTALL_PATH = /usr/local/cuda
SDK_INSTALL_PATH = /usr/local/cuda_sdk/C
VTK_INCLUDE_PATH = /usr/include/vtk-5.4

DEVICE  = cpu
#DEVICE  = gpu

#EXPAND  = Cartesian
EXPAND  = Spherical

CXX     = mpicxx -mpreferred-stack-boundary=4 -ggdb3 -Wall -Wextra -Winit-self -Wshadow -O3 -fPIC -fopenmp\
	-ffast-math -funroll-loops -fforce-addr -rdynamic -D_FILE_OFFSET_BITS=64\
	-I../include -I$(VTK_INCLUDE_PATH)
#CXX     = mpicxx -O2 -fPIC -openmp -I../include -I$(VTK_INCLUDE_PATH)
NVCC    = nvcc -Xcompiler -fopenmp --ptxas-options=-v -O3 -use_fast_math -arch=sm_13\
	-I../include -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc
LFLAGS  = -D$(DEVICE) -D$(EXPAND)
ifeq ($(DEVICE),gpu)
LFLAGS  += -L$(CUDA_INSTALL_PATH)/lib64 -L$(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lcutil_x86_64 -lstdc++ -ldl -lm
endif
#VFLAGS  = -lvtkCommon -lvtkCharts -DVTK
OBJECT  = ../kernel/$(DEVICE)$(EXPAND)Laplace.o ../kernel/$(DEVICE)$(EXPAND)BiotSavart.o\
	../kernel/$(DEVICE)$(EXPAND)Stretching.o ../kernel/$(DEVICE)$(EXPAND)Gaussian.o\
	../kernel/$(DEVICE)$(EXPAND)CoulombVdW.o

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)
.cu.o   :
	$(NVCC) -c $? -o $@ $(LFLAGS)
cleanall:
	rm -rf `find .. -name "*.o" -o -name "*.out*" -o -name "*.dat" -o -name "*.a" -o -name "*.sum"`
	rm -f ../unit_test/direct0*
commit  :
	hg commit
	hg push
	hg pull -u
save    :
	make cleanall
	tar zcvf ../../exafmm.tgz ../../exafmm
docs:
	doxygen Doxyfile
	cd docs/html ;tar czf ../../docs.tar *
	scp docs.tar $(EXAFMM_DOCS_USER)@barbagroup.bu.edu:~/
	ssh $(EXAFMM_DOCS_USER)@barbagroup.bu.edu 'tar -xmzf docs.tar -C /Library/WebServer/Documents/exafmm_docs/html/; rm docs.tar; chmod -R 775 /Library/WebServer/Documents/exafmm_docs/'
	rm -rf docs*

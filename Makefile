# These will move to coompilerOptions.py and some other tests
GNU_FLAGS    = -ggdb3 -Wall -Wextra -Winit-self -Wshadow -O3 -fPIC -fopenmp -ffast-math -funroll-loops -fforce-addr -rdynamic -D_FILE_OFFSET_BITS=64
INTEL_FLAGS  = -O2 -fPIC -openmp
NVIDIA_FLAGS = -Xcompiler -fopenmp --ptxas-options=-v -O3 -use_fast_math -arch=sm_21 -I$(CUDA_INSTALL_PATH)/include -I$(SDK_INSTALL_PATH)/common/inc
CUDA_LFLAGS  = -L$(CUDA_INSTALL_PATH)/lib64 -L$(SDK_INSTALL_PATH)/lib -lcuda -lcudart -lcutil_x86_64 -lstdc++ -ldl -lm
VTK_FLAGS    = -I/usr/include/vtk-5.6 -lvtkRendering -lvtkGraphics -lvtkFiltering -lvtkViews -lvtkCommon -lvtkWidgets -lvtkIO -DVTK

all:
	${MAKE} -C kernel lib

commit:
	hg commit
	hg push
	hg pull -u

save:
	hg st -un0 | xargs -0 rm
	tar zcf ../../exafmm.tgz ../../exafmm

docs:
	doxygen Doxyfile
	cd docs/html; tar czf ../../docs.tar *
	scp docs.tar pl:~/
	ssh pl 'tar -xmzf docs.tar -C /Library/WebServer/Documents/exafmm_docs/html/; rm docs.tar; chmod -R 775 /Library/WebServer/Documents/exafmm_docs/'
	rm -rf docs*

.PHONY: commit save docs

#
# Memo
#
#
# 080609 MD_PERIODIC_FIXED=1 seems stragne
#          ./emutest a 0  generates about 10% of error
#            this might be caused by insufficient dynamic range.
#          emulator does not work well
#
# 080609d ./emutest x 100 ---- seems ok
#         ./emutest x 10  ---- average error is 10^-3, which is large. Why?
#
# 080612  MD_SORT_ATYPEI does't work with AMBER. Last 128 atoms (?) seems strange.
#
# 080613  overlap with pthread in MR3calccoulomb_vdw_ij_exlist does't work. 
#         and it seems slow

# 100216  cutoff length can be specified by mr3_set_rcut2_ or vg_set_rcut2.
#         this is checked with c36a3xq2

.SUFFIXES: .cu .c .o

CC = gcc
C++ = g++
#CC = icc -openmp -xW -tpp7
#C++ = icpc -openmp -xW -tpp7
F90 = gfortran
F77 = gfortran
#F90 = ifort -xW -tpp7
#F77 = $(F90)

CUDA_SDK2=-DCUDA_SDK_2
#CUDA_SDK2=

CFLAGSCORE = -O $(FLAGS)
#CFLAGSCORE = -O0
#CFLAGSCORE += -g 
#CFLAGSCORE += -pg 
CFLAGSCORE += -pg -g
CFLAGSCORE += -ffast-math
#CFLAGSCORE += -fopenmp
#CFLAGSCORE += -ftree-vectorize -ftree-vectorizer-verbose=2 -msse3

INCLUDEDIR = ../../include
INCLUDEDEP1 = $(INCLUDEDIR)/mdgrape3.h
#CFLAGS = $(COMMONFLAGS)  $(CFLAGSCORE) -I$(INCLUDEDIR)
CFLAGS = $(CFLAGSCORE)
CFLAGS += $(CUDA_SDK2)
FFLAGS = $(CFLAGSCORE)
G80LIB = libmd3withg80.a
G80LIBVG = libmd3withg80vg.a
#TARGET = emutestvg $(G80LIB) $(G80LIBVG) sample_md3 forttest sample_md3vg forttestvg emutest
#TARGET = emutestvg sample_md3vg sample_md3
#TARGET = emutestvg sample_md3vg r1.g80emu rsqrt.g80emu
TARGET = emutestvg sample_md3vg
#TARGET = sample_md3
#GPUCOREOBJ = gpuvdwforce_kernel.o gpuvdwforce_mother.o \
#	gpuvdwpot_kernel.o gpuvdwpot_mother.o \
#	gpucoulombforce_kernel.o gpucoulombforce_mother.o \
#	gpucoulombpot_kernel.o gpucoulombpot_mother.o \
#	gpuewreforce_kernel.o gpuewreforce_mother.o \
#	gpuewrepot_kernel.o gpuewrepot_mother.o \
#	gpuewwvforce_kernel.o gpuewwvforce_mother.o \
#	gpuewwvpot_kernel.o gpuewwvpot_mother.o 
GPUCOREOBJ = gpuvdwforce_mother.o gpuvdwpot_mother.o gpucoulombforce_mother.o \
	gpucoulombpot_mother.o gpuewreforce_mother.o gpuewrepot_mother.o \
	gpuewwvforce_mother.o gpuewwvpot_mother.o 

CUDA_INSTALL_PATH=/usr/local/cuda
NVCC       = $(CUDA_INSTALL_PATH)/bin/nvcc 
NVCC_FLAGS = -cuda --host-compilation C 
#NVCC_FLAGS = --compile
NVCC_FLAGS += -use_fast_math
NVCC_FLAGS += --ptxas-options=-v
#NVCC_FLAGS += -deviceemu
#NVCC_FLAGS += -gpu
#NVCC_FLAGS += -ptx
#NVCC_FLAGS += -cubin
NVCC_FLAGS += -arch sm_11
#NVCC_FLAGS += -arch sm_13
#NVCC_FLAGS += -arch sm_20
NVCC_FLAGS += $(CUDA_SDK2)
NVCC_FLAGS += -O
# for debug with SDK
#CUDAINCLUDES  = -I. -I$(CUDA_INSTALL_PATH)/include -I$(CUDA_SDK_PATH)/common/inc
CUDAINCLUDES  = -I. -I$(CUDA_INSTALL_PATH)/include
#CUDALIB = $(CUDA_SDK_PATH)/lib/libcutil.a \
#		-L$(CUDA_INSTALL_PATH)/lib -lcuda -lcudart -lGLU -lGL -lm
#CUDALIB = $(CUDA_SDK_PATH)/lib/libcutil.a \
#		-L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lGLU -lGL -lm
# for debug with SDK
#CUDALIB = $(CUDA_SDK_PATH)/lib/libcutil_x86_64.a \
#		-L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lGLU -lGL -lm
CUDALIB = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lGLU -lGL -lm -lpthread

.cu.c:
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
.c.o:
	$(C++) -c $(CFLAGS) $<
#	$(CC) -c $(CFLAGS) $<

all : $(TARGET)

#gpucoulombforce_mother.o : gpucoulombforce_mother.ptx

$(G80LIB) : mr3.o $(GPUCOREOBJ)
	ar ruv $@ $?
	ranlib $@
mr3.o : mr3.c mr3_host.c
	$(CC) -c $(CFLAGS) $<
vtgrape.o : vtgrape.c vtgrape.h md.h md.c mr3_host.c \
            vtgrapeproto.h vtgrape_mr3.c vtgrape_cellindex.c \
            vtgrape_dummy.c sock_lib.c vtgrape_mixed.c mddriver.c 
	$(CC) -c $(CFLAGS) $<
#vtgrape.o : vtgrape.c vtgrape.h md.h md.c mr3_host.c md_mother.o
#	$(CC) -c $(CFLAGS) vtgrape.c md_mother.o

#func.o : func.c++ md.h
#	$(C++) -c $(CFLAGS) $<

forttest : forttest.f90 $(G80LIB)
	$(F90) -o $@ $< $(G80LIB) $(CUDALIB)
sample_md3 : sample_md3.f $(G80LIB)
	$(F77) $(FFLAGS) -o $@ $< $(G80LIB) $(CUDALIB)
forttestvg : forttest.f90 vtgrape.o md_mother.o
	$(F90) -o $@ $< vtgrape.o $(CUDALIB) md_mother.o
sample_md3vg : sample_md3.f $(G80LIBVG)
	$(F77) $(FFLAGS) -o $@ $< $(G80LIBVG) $(CUDALIB) 

vtgrapeproto.h : vtgrape.c
	cproto vtgrape.c >> vtgrapeproto.h

emutest : emutest.c $(G80LIB)
	$(CC) $(CFLAGS) -o $@ $< $(G80LIB) $(CUDALIB) 
emutestvg : emutest.c $(G80LIBVG)
	$(CC) $(CFLAGS) -o $@ $< $(G80LIBVG) $(CUDALIB) 

r1.g80emu : emutestvg initialize.sh
	./initialize.sh

#md_mother.o : md_mother.cu md_kernel.cu
#	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@

$(G80LIBVG) : md_mother.o vtgrape.o 
	ar ruv $@ $?
	ranlib $@

gpucoulombforce_mother.o : gpucoulombforce_mother.cu gpucoulombforce_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpucoulombpot_mother.o : gpucoulombpot_mother.cu gpucoulombpot_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpuvdwforce_mother.o : gpuvdwforce_mother.cu gpuvdwforce_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpuvdwpot_mother.o : gpuvdwpot_mother.cu gpuvdwpot_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpuewreforce_mother.o : gpuewreforce_mother.cu gpuewreforce_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpuewrepot_mother.o : gpuewrepot_mother.cu gpucoulombpot_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpuewwvforce_mother.o : gpuewwvforce_mother.cu gpuewwvforce_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@
gpuewwvpot_mother.o : gpuewwvpot_mother.cu gpuewwvpot_kernel.cu
	$(NVCC) $(NVCC_FLAGS) $(CUDAINCLUDES) $< -o $@

dist :
	tar chvfz md3withg80.tar.gz Makefile Makefile_intel *.[ch] *.cu *.f *.sh *.md5sum_org *.c++

clean : 
	rm -f $(TARGET) *.o *~ kernel.c gpu*.c md_mother.c \
		$(G80LIB) $(G80LIBVG) *.md5sum

# memo
#
# 070928 Ewald wave potential is calculated only total potential
#        Ewald wave force accuracy is only 10^-2.2
#        vdw force and pot accuracy is not good
# 080924 potential result in MR3calccoulomb_nlist_ij_emu and
#        MR3calcvdw_nlist_ij_emu are only force[i*3] not 
#        force[i*3+1] nor force[i*3+2].

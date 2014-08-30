.SUFFIXES: .cxx .o

CXX	= mpicxx -ggdb3 -Wall -O3 -fopenmp
#CXX	= CC -xHOST -O3 -openmp -funroll-loops
#CXX	= mpixlcxx_r -qarch=450 -qtune=450 -O3
#CXX	= mpiFCCpx -Kfast,openmp

FC 	= mpif90 -ggdb3 -O3 -msse4a -ffast-math -funroll-loops -cpp
LFLAGS 	=

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)

serial	: main.cxx
	$(CXX) $? $(LFLAGS) -DDO_P2P -DSerial
	./a.out

parallel: main.cxx
	$(CXX) $? $(LFLAGS) -DDO_P2P #-DPRINT_COMM
	mpirun -np 2 ./a.out

wrapper:
	make libfmm.a
	$(FC) test_wrapper.f90 -ggdb3 -Wall -fbounds-check -L. -lfmm -ltbb -lstdc++ -lgomp
	mpirun -np 1 ./a.out ../wrappers/ala3.fmm
	make cleandat
	mpirun -np 1 ./a.out ../wrappers/water.fmm 10

libfmm.a: wrapper.o
	ar ruv libfmm.a $?
	ranlib libfmm.a

cleanjob:
	find . -name "job.sh*" | xargs rm -rf
clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
cleandat:
	find . -name "*.dat" | xargs rm -rf
commit  :
	hg commit
	hg push
	hg pull -u
save    :
	make clean
	cd .. && tar zcvf fmm.tgz fmm_cppmd
revert  :
	hg revert --all
	rm -rf `find . -name "*.orig"`

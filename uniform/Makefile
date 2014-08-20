CXX	= mpicxx -ggdb3 -Wall -O3 -fopenmp
#CXX	= CC -xHOST -O3 -openmp -funroll-loops
#CXX	= mpixlcxx_r -qarch=450 -qtune=450 -O3
#CXX	= mpiFCCpx -Kfast,openmp

LFLAGS =

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)

serial	: main.cxx
	$(CXX) $? $(LFLAGS) -DDO_P2P -DSerial
	./a.out

parallel: main.cxx
	$(CXX) $? $(LFLAGS) -DDO_P2P #-DPRINT_COMM
	mpirun -np 2 ./a.out

cleanjob:
	find . -name "job.sh*" | xargs rm -rf	
clean:
	find . -name "*.o" -o -name "*.out*" -o -name "*.dat" | xargs rm -rf
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

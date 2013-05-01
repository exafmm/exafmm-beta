include ../Makefile.include

MR3 = mr3_host.o mr3.o

libserial_coulombVdW.a: serial_coulombVdW.o $(OBJECT)
	ar ruv libserial_coulombVdW.a $?
	ranlib libserial_coulombVdW.a

libserial_Ewald.a: serial_Ewald.o $(OBJECT)
	ar ruv libserial_Ewald.a $?
	ranlib libserial_Ewald.a
	mv libserial_Ewald.a libserial_coulombVdW.a

libparallel_coulombVdW.a: parallel_coulombVdW.o $(OBJECT)
	ar ruv libparallel_coulombVdW.a $?
	ranlib libparallel_coulombVdW.a

test_serial_coulombVdW: test_serial_coulombVdW.cxx
	rm -f libserial_coulombVdW.a
	make libserial_coulombVdW.a
	$(CXX) $? -L. -lserial_coulombVdW $(LFLAGS)
	./a.out

test_parallel_coulombVdW: test_parallel_coulombVdW.cxx
	rm -f libparallel_coulombVdW.a
	make libparallel_coulombVdW.a
	$(CXX) $? -L. -lparallel_coulombVdW -lmd3withg80vg $(LFLAGS)
	mpirun -np 2 ./a.out

ewald: ewald.cxx $(MR3)
	rm -f libserial_coulombVdW.a
	make libserial_coulombVdW.a
	$(CXX) $? -L. -lserial_coulombVdW $(LFLAGS)
	./a.out 0 12
	./a.out 1 12
	./a.out 2 12
	./a.out 3 12

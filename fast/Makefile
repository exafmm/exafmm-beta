include ../Makefile

a	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT
	./a.out

b	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT
	./a.out

d	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out


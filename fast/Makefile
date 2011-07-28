include ../Makefile

a	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT
	./a.out

a2	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DCART
	./a.out

a3	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DSPHE
	./a.out

b	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

b2	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DCART
	./a.out

b3	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DSPHE
	./a.out

c	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT
	./a.out

c2	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DCART
	./a.out

c3	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DSPHE
	./a.out

d	: main.cc
	make cleanall
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out


include ../Makefile

a	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DCART -DIEQJ
	./a.out

a2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DSPHE -DIEQJ
	./a.out

b	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DCART -DIEQJ
	./a.out

b2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DSPHE -DIEQJ
	./a.out

c	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DCART -DIEQJ
	./a.out

c2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DSPHE -DIEQJ
	./a.out

d	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

f	: fdgl.cxx
	make cleanall
	$(CXX) $? -lvtkCharts
	./a.out

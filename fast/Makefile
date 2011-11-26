include ../Makefile

a	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DCART
	./a.out

a2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DSPHE
	./a.out

b	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DCART
	./a.out

b2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DSPHE
	./a.out

c	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DCART
	./a.out

c2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DSPHE
	./a.out

d	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DBUILD -DCART
	./a.out

f	: fdgl.cxx
	make cleanall
	$(CXX) $? -lvtkCommon -lvtkCharts -lvtkFiltering -lvtkViews -lvtkRendering -DMANY -DCART
	./a.out

include ../Makefile

a	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DCartesian
	./a.out

a2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT -DSpherical
	./a.out

b	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DCartesian
	./a.out

b2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DSpherical
	./a.out

c	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DCartesian
	./a.out

c2	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT -DSpherical
	./a.out

d	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DBUILD -DCartesian
	./a.out

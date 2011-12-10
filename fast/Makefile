include ../Makefile
#LFLAGS	+= -DIneJ
#LFLAGS	+= -DBOTTOMUP

a	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DDIRECT
	./a.out

b	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DMANY -DDIRECT
	./a.out

d	: main.cxx
	make cleanall
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

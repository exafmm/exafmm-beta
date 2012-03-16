include ../Makefile.include
LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE

a	: main.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	./a.out

b	: main.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

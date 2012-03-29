include ../Makefile.include
#LFLAGS	+= -DIneJ
#LFLAGS	+= -DBOTTOMUP
#LFLAGS	+= -DHYBRID
#LFLAGS	+= -DTREECODE
#LFLAGS += -DRECURSIVE -DCOMMAND_LINE

a.out : main.cxx
	$(CXX) -o $@ $< $(LFLAGS)

a	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS)
	./a.out

b	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DMANY
	./a.out

c	: main.cxx
	make clean
	$(CXX) $? $(LFLAGS) -DBUILD
	./a.out

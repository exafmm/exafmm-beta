#CXX	= mpicxx -ggdb3 -Wall -Wextra -O3 -ffast-math -funroll-loops -fopenmp
CXX    = icpc -Wall -xHOST -O3 -funroll-loops -finline-functions -ansi-alias -openmp -I/usr/include/x86_64-linux-gnu/c++/4.8 # Intel
MAKEFLAGS += --no-print-directory

LFLAGS += -lsakura
LFLAGS += -lstencil

all	:
	make clean
	make compile
	make run

compile	: main.cxx
	$(CXX) $? $(LFLAGS)

run	:
	./a.out

clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf

CXX	= mpicxx -ggdb3 -Wall -Wextra -O3 -ffast-math -funroll-loops -fopenmp
MAKEFLAGS += --no-print-directory

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

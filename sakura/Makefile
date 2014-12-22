#CXX	= mpicxx -ggdb3 -Wall -Wextra -O3 -ffast-math -funroll-loops -fopenmp
CXX    = icpc -Wall -xHOST -O3 -funroll-loops -finline-functions -ansi-alias -openmp -std=gnu99 -I/usr/include/x86_64-linux-gnu/c++/4.8 # Intel
MAKEFLAGS += --no-print-directory

LFLAGS	+=

all	:
	make clean
	make compile
	make run

compile	: test_sakura.c datasets.c quantization.c build_tree_cilk.c data_rearrange.c morton_encoding.c cilk_linked_list.c verification.c wrapper.c utils.c
	$(CXX) $? $(LFLAGS)

run	:
	./a.out 10000000 2 64 1 8

clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf

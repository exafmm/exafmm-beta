.SUFFIXES: .c .o

CXX	= icpc -Wall -xHOST -O3 -funroll-loops -finline-functions -ansi-alias -openmp -std=gnu99 -I/usr/include/x86_64-linux-gnu/c++/4.8
SOURCES	= test_sakura.c datasets.c quantization.c build_tree_cilk.c data_rearrange.c morton_encoding.c cilk_linked_list.c verification.c wrapper.c utils.c
OBJECTS = $(SOURCES:.c=.o)
LFLAGS	+=
MAKEFLAGS += --no-print-directory

.c.o	:
	$(CXX) -c $^ -o $@

all	:
	make compile
	make run

compile	: $(OBJECTS)
	$(CXX) $? $(LFLAGS)

run	:
	./a.out 10000000 2 64 1 8

clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf

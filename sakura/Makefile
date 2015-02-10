.SUFFIXES: .c .o

CXX	= icc -Wall -xHOST -O3 -funroll-loops -finline-functions -ansi-alias -openmp -std=gnu99 -I/usr/include/x86_64-linux-gnu/c++/4.8
SOURCES	= test_sakura.c datasets.c encode.c sort.c build_tree.c linked_list.c kernels.c
OBJECTS = $(SOURCES:.c=.o)
LFLAGS	+= -DTEST
MAKEFLAGS += --no-print-directory

.c.o	:
	$(CXX) -c $^ -o $@ $(LFLAGS)

all	:
	make clean
	make compile
	make run

compile	: $(OBJECTS)
	$(CXX) $? $(LFLAGS)

run	:
	./a.out 4 2 1 1
#	./a.out 100000 1 64 4

clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf

revert:
	hg revert --all
	find . -name "*.orig" | xargs rm -f

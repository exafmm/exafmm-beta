.SUFFIXES: .c .o

CXX	= icc -Wall -xHOST -O1 -funroll-loops -finline-functions -ansi-alias -openmp -std=gnu99 -I/usr/include/x86_64-linux-gnu/c++/4.8
SOURCES	= test_sakura.c datasets.c encode.c sort.c build_tree.c linked_list.c kernels.c
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
	./a.out 1000000 2 64 4

clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf

revert:
	hg revert --all
	find . -name "*.orig" | xargs rm -f

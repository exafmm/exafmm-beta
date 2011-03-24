include ../Makefile

static: $(KERNEL).o $(OBJECT)
	ar ruv ../lib/lib$(DEVICE)$(KERNEL).a $?
	ranlib ../lib/lib$(DEVICE)$(KERNEL).a

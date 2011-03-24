include ../Makefile

static: $(KERNEL).o $(OBJECT)
	ar ruv ../lib/lib$(DEVICE)$(KERNEL).a $?
	ranlib ../lib/lib$(DEVICE)$(KERNEL).a
dynamic: $(KERNEL).o $(OBJECT)
	$(CXX) $^ -shared -o ../lib/lib$(DEVICE)$(KERNEL).so
dynamic2: $(KERNEL).o $(OBJECT)
	ld -r $(KERNEL).o $(OBJECT) -o $(DEVICE)$(KERNEL).o
	$(CXX) $(DEVICE)$(KERNEL).o -shared -o ../lib/lib$(DEVICE)$(KERNEL).so

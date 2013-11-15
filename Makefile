help:
	@echo "Please read the README file in the root directory."
clean:
	find . -name "*.o" -o -name "*.out*" -o -name "*.mod" | xargs rm -rf
cleandat:
	find . -name "*.dat" -o -name "*.pdb" | xargs rm -f
cleanlib:
	find . -name "*.a" -o -name "*.so" | xargs rm -f
cleanall:
	make clean
	make cleandat
	make cleanlib
commit  :
	hg commit
	hg push
	hg pull -u
save    :
	make cleanall
	cd .. && tar zcvf exafmm-dev.tgz exafmm-dev
revert	:
	hg revert --all
	find . -name "*.orig" | xargs rm -f

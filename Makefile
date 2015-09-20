help:
	@echo "Please read the README file in the root directory."
buildbot:
	buildbot try --connect=pb --master=localhost:9989 --username=exafmm --passwd=exafmm --vc=git
clean:
	@find . -name "*.o" -o -name "*.out*" -o -name "*.mod" | xargs rm -rf
cleandat:
	@find . -name "*.dat" -o -name "*.pdb" -o -name "*_restart" | xargs rm -f
cleanlib:
	@find . -name "*.a" -o -name "*.so" | xargs rm -f
cleanall:
	make clean
	make cleandat
	make cleanlib
revert:
	git reset --hard HEAD
save:
	make cleanall
	cd .. && tar zcvf exafmm.tgz exafmm
tags:
	ctags -eR .

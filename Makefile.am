SUBDIRS = examples wrappers helmholtz

help:
	@echo "Please read the README file in the root directory."
buildbot:
	buildbot try --connect=pb --master=131.112.5.69:443 --username=exafmm --passwd=exafmm --vc=git
cleanbin:
	@find . -name "*.o" -o -name "*.out*" -o -name "*.mod" | xargs rm -rf
cleandat:
	@find . -name "*.dat" -o -name "*.pdb" -o -name "*.svg" -o -name "*_restart" | xargs rm -f
cleanlib:
	@find . -name "*.a" -o -name "*.so" | xargs rm -f
cleanall:
	make cleanbin
	make cleandat
	make cleanlib
revert:
	git reset --hard HEAD
save:
	make cleanall
	cd .. && tar zcvf exafmm.tgz exafmm

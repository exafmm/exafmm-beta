help:
	@echo "Help documentation will be available soon...\n"
clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
cleandat:
	find . -name "*.dat" -o -name "*.dot" -o -name "*.svg" | xargs rm -rf
cleanlib:
	find . -name "*.a" -o -name "*.so" | xargs rm -rf
cleanall:
	make clean
	make cleandat
commit  :
	hg commit
	hg push
	hg pull -u
save    :
	make cleanall
	cd .. && tar zcvf exafmm.tgz exafmm
revert	:
	hg revert --all
	rm -rf `find . -name "*.orig"`

help:
	@echo "Help documentation will be available soon.\n"
clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
cleandat:
	find . -name "*.dat" -o -name "*.dot" -o -name "*.svg" | xargs rm -rf
cleanlib:
	find . -name "*.a" -o -name "*.so" | xargs rm -rf
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
	cd .. && tar zcvf exafmm.tgz exafmm
revert	:
	hg revert --all
	rm -rf `find . -name "*.orig"`
docs:
	doxygen Doxyfile
	cd docs/html; tar zcf ../../docs.tgz *
	scp docs.tgz pl:
	ssh pl 'tar -zxf docs.tgz -C /Library/WebServer/Documents/exafmm_docs/html/; rm docs.tgz; chmod -R 775 /Library/WebServer/Documents/exafmm_docs/'
	rm -rf docs*

cleanall:
	find . -name "*.o" -o -name "*.out*" -o -name "*.dat" -o -name "*.dot" -o -name "*.svg" | xargs rm -rf
commit  :
	hg commit
	hg push
	hg pull -u
save    :
	make cleanall
	cd .. && tar zcvf exafmm.tgz exafmm
revert	:
	hg revert --all
	rm -rf `find -name "*.orig"`
docs:
	doxygen Doxyfile
	cd docs/html; tar czf ../../docs.tar *
	scp docs.tar pl:~/
	ssh pl 'tar -xmzf docs.tar -C /Library/WebServer/Documents/exafmm_docs/html/; rm docs.tar; chmod -R 775 /Library/WebServer/Documents/exafmm_docs/'
	rm -rf docs*

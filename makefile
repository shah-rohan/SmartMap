default :
	$(MAKE) -C gzstream default
	$(MAKE) -C Default all

clean : 
	$(MAKE) -C gzstream cleanmore
	$(MAKE) -C Default clean

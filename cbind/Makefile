
include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)
LIBNAME=libmld_cbind.a	

all: lib
lib: mlprecd
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *.h $(INCDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)	


mlprecd:
	$(MAKE) -C mlprec lib LIBNAME=$(LIBNAME) 


clean: 
	$(MAKE) -C mlprec clean


veryclean: clean
	cd test/pargen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h

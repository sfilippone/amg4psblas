
include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)
LIBNAME=libamg_cbind.a	

objs: amgprecd
lib: objs
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *.h $(INCDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)


amgprecd:
	cd amgprec && $(MAKE)  objs LIBNAME=$(LIBNAME)


clean:
	cd amgprec &&$(MAKE)  clean


veryclean: clean
	cd test/pargen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h

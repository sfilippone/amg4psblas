
include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)
LIBNAME=libamg_cbind.a	
AMGFDEFINES=-DPSB_HAVE_LAPACK -DPSB_HAVE_FLUSH_STMT -DPSB_MPI_MOD  $(PSBFDEFINES) $(FCUDEFINES)
FDEFINES=$(AMGFDEFINES) 

objs: amgprecd

lib: objs
	cd amgprec && $(MAKE)  lib LIBNAME=$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)


amgprecd:
	cd amgprec && $(MAKE)  objs


clean:
	cd amgprec &&$(MAKE)  clean


veryclean: clean
	cd test/pargen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h


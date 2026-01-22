include ../../../Make.inc
#
# Libraries used
#
LIBDIR=../../../lib/
INCLUDEDIR=../../../include/
HERE=.

FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBLAS_INCDIR)
#PSBLAS_LIBS= -L$(PSBLAS_LIBDIR) -L$(LIBDIR)  $(CPSBLAS_LIB) $(PSBLAS_LIB)
# -lpsb_linsolve_cbind -lpsb_prec_cbind -lpsb_base_cbind 
PSBC_LIBS= -L$(PSBLAS_LIBDIR) -lpsb_cbind  -lpsb_linsolve -lpsb_prec 
AMGC_LIBS=-L$(LIBDIR) -lamg_cbind  -lamg_prec 
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).
CINCLUDES=-I$(LIBDIR) -I$(PSBLAS_INCDIR) -I$(INCLUDEDIR)

EXEDIR=./runs
# 20110404 specifying UMFLIBS here is not portable!
#UMFLIBS=-lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd -L/usr/include/suitesparse
#UMFFLAGS=-DHave_UMF_ -I/usr/include/suitesparse

 all:  amgec amgecgpu

amgec: amgec.o
	$(MPFC) amgec.o  -o amgec $(AMGC_LIBS)  $(PSBC_LIBS) $(PSBCLDLIBS) $(PSBLAS_LIBS) \
	  $(UMFLIBS) $(PSBLDLIBS) $(AMGLDLIBS) $(PSBGPULDLIBS) $(LDLIBS) -lm -lgfortran -fopenmp
# \
#	-lifcore -lifcoremt -lguide -limf -lirc -lintlc -lcxaguard -L/opt/intel/fc/10.0.023/lib/ -lm

	/bin/mv amgec $(EXEDIR)

amgecgpu: amgecgpu.o
	$(MPFC) amgecgpu.o  -o amgecgpu $(AMGC_LIBS)  $(PSBC_LIBS) $(PSBCLDLIBS) $(PSBLAS_LIBS) \
	  $(UMFLIBS) $(PSBLDLIBS) $(AMGLDLIBS) $(PSBGPULDLIBS) $(LDLIBS) -lm -lgfortran -fopenmp
	/bin/mv amgecgpu $(EXEDIR)

.f90.o:
	$(MPFC) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<
.c.o:
	$(MPCC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<


clean: 
	/bin/rm -f amgec.o  $(EXEDIR)/amgec $(EXEDIR)/amgecgpu
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)

tests: all
	cd runs ; ./amgec < amge.inp
	cd runs ; ./amgecgpu < amgegpu.inp



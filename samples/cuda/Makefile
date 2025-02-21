AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
# To compile these example you need to have also installed PSBLAS-EXT, SPGPU and
# CUDA. Please refer to the documentation of AMG4PSBLAS and PSBLAX-EXT for the
# full details.
# This define the variables : PSBGPULDLIBS, LIBRSB_INCLUDES, and CUDA_INCLUDES
#PSBEXT_DIR= #< FILL-ME >
#include $(PSBEXT_DIR)/include/Make.inc.ext
include $(PSBLASDIR)/include/Make.inc.psblas
PSBLIBDIR=$(PSBLASDIR)/lib/
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) \
$(PSBLAS_INCLUDES) $(FIFLAG). 
AMG_LIB=-L$(AMGLIBDIR) -lamg_prec 
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_linsolve  -lpsb_prec -lpsb_util -lpsb_ext  \
	 -lpsb_base 

LINKOPT=-lstdc++

DGOBJS=amg_dexample_cuda.o data_input.o

EXEDIR=./runs

all: amg_dexample_cuda


amg_dexample_cuda: $(DGOBJS)
	$(FLINK) $(LINKOPT) $(DGOBJS) -o amg_dexample_cuda \
	$(AMG_LIB) $(PSBGPULDLIBS) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv amg_dexample_cuda $(EXEDIR)

amg_dexample_cuda.o: data_input.o 

clean:
	/bin/rm -f $(DGOBJS) *$(.mod)  $(EXEDIR)/amg_dexample_cuda

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

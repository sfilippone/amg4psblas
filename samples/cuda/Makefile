AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
# This code will include the appropriate defines and libraries
# if you have compiled the main PSBLAS library with CUDA support
include $(PSBLASDIR)/include/Make.inc.psblas
PSBLIBDIR=$(PSBLASDIR)/lib/
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) \
$(PSBLAS_INCLUDES) $(FIFLAG). 
AMG_LIB=-L$(AMGLIBDIR) -lamg_prec 
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_linsolve  -lpsb_prec -lpsb_util -lpsb_ext  \
	 -lpsb_base 
FDEFINES=$(AMGFDEFINES) $(FCUDEFINES)

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

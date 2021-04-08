AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).
LINKOPT=

DGOBJS=amg_dexample_gpu.o data_input.o amg_dpde_mod.o

EXEDIR=./runs

all: amg_dexample_gpu


amg_dexample_gpu: $(DGOBJS)
	$(FLINK) $(LINKOPT) $(DGOBJS) -o amg_dexample_gpu \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_dexample_gpu $(EXEDIR) 

amg_dexample_gpu.o: data_input.o amg_dpde_mod.o

clean: 
	/bin/rm -f $(DGOBJS) *$(.mod)  $(EXEDIR)/amg_dexample_gpu

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


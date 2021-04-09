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
PSBEXT_DIR= #< FILL-ME >
include $(PSBEXT_DIR)/include/Make.inc.ext

FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) \
$(PSBLAS_INCLUDES) $(FIFLAG). $(FMFLAG)$(SPGPU_INCDIR) \
$(FMFLAG)$(LIBRSB_INCLUDES) $(FMFLAG)$(CUDA_INCLUDES) \

LINKOPT=$(CINCLUDES) $(CDEFINES) $(CUDEFINES) 

DGOBJS=amg_dexample_gpu.o data_input.o amg_dpde_mod.o

EXEDIR=./runs

all: amg_dexample_gpu


amg_dexample_gpu: $(DGOBJS)
	$(FLINK) $(LINKOPT) $(DGOBJS) -o amg_dexample_gpu \
	$(AMG_LIBS) $(PSBLAS_LIBS) $(PSBGPULDLIBS) $(LDLIBS)
	/bin/mv amg_dexample_gpu $(EXEDIR)

amg_dexample_gpu.o: data_input.o amg_dpde_mod.o

clean:
	/bin/rm -f $(DGOBJS) *$(.mod)  $(EXEDIR)/amg_dexample_gpu

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

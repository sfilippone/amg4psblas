AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).
LINKOPT=

DMOBJS=amg_dexample_ml.o data_input.o amg_dpde_mod.o
D1OBJS=amg_dexample_1lev.o data_input.o amg_dpde_mod.o
SMOBJS=amg_sexample_ml.o data_input.o amg_spde_mod.o
S1OBJS=amg_sexample_1lev.o data_input.o amg_spde_mod.o

EXEDIR=./runs

all: amg_sexample_ml amg_sexample_1lev amg_dexample_ml amg_dexample_1lev


amg_dexample_ml: $(DMOBJS)
	$(FLINK) $(LINKOPT) $(DMOBJS) -o amg_dexample_ml \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_dexample_ml $(EXEDIR) 

amg_dexample_1lev: $(D1OBJS)
	$(FLINK) $(LINKOPT) $(D1OBJS) -o amg_dexample_1lev \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_dexample_1lev $(EXEDIR) 

amg_dexample_ml.o: data_input.o amg_dpde_mod.o
amg_dexample_1lev.o: data_input.o amg_dpde_mod.o

amg_sexample_ml: $(SMOBJS)
	$(FLINK) $(LINKOPT) $(SMOBJS) -o amg_sexample_ml \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_sexample_ml $(EXEDIR) 

amg_sexample_1lev: $(S1OBJS)
	$(FLINK) $(LINKOPT) $(S1OBJS) -o amg_sexample_1lev \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_sexample_1lev $(EXEDIR) 

amg_sexample_ml.o: data_input.o  amg_spde_mod.o
amg_sexample_1lev.o: data_input.o amg_spde_mod.o

clean: 
	/bin/rm -f $(DMOBJS) $(D1OBJS) $(SMOBJS) $(S1OBJS) \
        *$(.mod)  $(EXEDIR)/amg_dexample_ml $(EXEDIR)/amg_dexample_1lev\
	 $(EXEDIR)/amg_sexample_ml $(EXEDIR)/amg_sexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


AMGDIR=../../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).
LINKOPT=

DMOBJS=amg_dexample_ml.o data_input.o
D1OBJS=amg_dexample_1lev.o data_input.o
ZMOBJS=amg_zexample_ml.o data_input.o
Z1OBJS=amg_zexample_1lev.o data_input.o
SMOBJS=amg_sexample_ml.o data_input.o
S1OBJS=amg_sexample_1lev.o data_input.o
CMOBJS=amg_cexample_ml.o data_input.o
C1OBJS=amg_cexample_1lev.o data_input.o

EXEDIR=./runs

all: amg_dexample_ml amg_dexample_1lev amg_zexample_ml amg_zexample_1lev\
	amg_sexample_ml amg_sexample_1lev amg_cexample_ml amg_cexample_1lev

amg_dexample_ml: $(DMOBJS)
	$(FLINK) $(LINKOPT) $(DMOBJS) -o amg_dexample_ml \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_dexample_ml $(EXEDIR) 

amg_dexample_1lev: $(D1OBJS)
	$(FLINK) $(LINKOPT) $(D1OBJS) -o amg_dexample_1lev \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_dexample_1lev $(EXEDIR) 

amg_dexample_ml.o: data_input.o
amg_dexample_1lev.o: data_input.o

amg_zexample_ml: $(ZMOBJS)
	$(FLINK) $(LINKOPT) $(ZMOBJS) -o amg_zexample_ml \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_zexample_ml $(EXEDIR) 

amg_zexample_1lev: $(Z1OBJS)
	$(FLINK) $(LINKOPT) $(Z1OBJS) -o amg_zexample_1lev \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_zexample_1lev $(EXEDIR) 

amg_zexample_ml.o: data_input.o
amg_zexample_1lev.o: data_input.o



amg_sexample_ml: $(SMOBJS)
	$(FLINK) $(LINKOPT) $(SMOBJS) -o amg_sexample_ml \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_sexample_ml $(EXEDIR) 

amg_sexample_1lev: $(S1OBJS)
	$(FLINK) $(LINKOPT) $(S1OBJS) -o amg_sexample_1lev \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_sexample_1lev $(EXEDIR) 

amg_sexample_ml.o: data_input.o
amg_sexample_1lev.o: data_input.o

amg_cexample_ml: $(CMOBJS)
	$(FLINK) $(LINKOPT) $(CMOBJS) -o amg_cexample_ml \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_cexample_ml $(EXEDIR) 

amg_cexample_1lev: $(C1OBJS)
	$(FLINK) $(LINKOPT) $(C1OBJS) -o amg_cexample_1lev \
	$(AMG_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv amg_cexample_1lev $(EXEDIR) 

amg_cexample_ml.o: data_input.o
amg_cexample_1lev.o: data_input.o


clean: 
	/bin/rm -f *$(.mod) \
	$(DMOBJS) $(D1OBJS) $(ZMOBJS) $(Z1OBJS) \
        $(EXEDIR)/amg_dexample_ml $(EXEDIR)/amg_dexample_1lev \
	$(EXEDIR)/amg_zexample_ml $(EXEDIR)/amg_zexample_1lev \
	$(SMOBJS) $(S1OBJS) $(CMOBJS) $(C1OBJS) \
        $(EXEDIR)/amg_sexample_ml $(EXEDIR)/amg_sexample_1lev \
	$(EXEDIR)/amg_cexample_ml $(EXEDIR)/amg_cexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


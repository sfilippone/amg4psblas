MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.amg4psblas
MLDMODDIR=$(MLDDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDMODDIR) $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).
LINKOPT=

DMOBJS=mld_dexample_ml.o data_input.o
D1OBJS=mld_dexample_1lev.o data_input.o
ZMOBJS=mld_zexample_ml.o data_input.o
Z1OBJS=mld_zexample_1lev.o data_input.o
SMOBJS=mld_sexample_ml.o data_input.o
S1OBJS=mld_sexample_1lev.o data_input.o
CMOBJS=mld_cexample_ml.o data_input.o
C1OBJS=mld_cexample_1lev.o data_input.o

EXEDIR=./runs

all: mld_dexample_ml mld_dexample_1lev mld_zexample_ml mld_zexample_1lev\
	mld_sexample_ml mld_sexample_1lev mld_cexample_ml mld_cexample_1lev

mld_dexample_ml: $(DMOBJS)
	$(FLINK) $(LINKOPT) $(DMOBJS) -o mld_dexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_dexample_ml $(EXEDIR) 

mld_dexample_1lev: $(D1OBJS)
	$(FLINK) $(LINKOPT) $(D1OBJS) -o mld_dexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_dexample_1lev $(EXEDIR) 

mld_dexample_ml.o: data_input.o
mld_dexample_1lev.o: data_input.o

mld_zexample_ml: $(ZMOBJS)
	$(FLINK) $(LINKOPT) $(ZMOBJS) -o mld_zexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_zexample_ml $(EXEDIR) 

mld_zexample_1lev: $(Z1OBJS)
	$(FLINK) $(LINKOPT) $(Z1OBJS) -o mld_zexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_zexample_1lev $(EXEDIR) 

mld_zexample_ml.o: data_input.o
mld_zexample_1lev.o: data_input.o



mld_sexample_ml: $(SMOBJS)
	$(FLINK) $(LINKOPT) $(SMOBJS) -o mld_sexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_sexample_ml $(EXEDIR) 

mld_sexample_1lev: $(S1OBJS)
	$(FLINK) $(LINKOPT) $(S1OBJS) -o mld_sexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_sexample_1lev $(EXEDIR) 

mld_sexample_ml.o: data_input.o
mld_sexample_1lev.o: data_input.o

mld_cexample_ml: $(CMOBJS)
	$(FLINK) $(LINKOPT) $(CMOBJS) -o mld_cexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_cexample_ml $(EXEDIR) 

mld_cexample_1lev: $(C1OBJS)
	$(FLINK) $(LINKOPT) $(C1OBJS) -o mld_cexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_cexample_1lev $(EXEDIR) 

mld_cexample_ml.o: data_input.o
mld_cexample_1lev.o: data_input.o


clean: 
	/bin/rm -f *$(.mod) \
	$(DMOBJS) $(D1OBJS) $(ZMOBJS) $(Z1OBJS) \
        $(EXEDIR)/mld_dexample_ml $(EXEDIR)/mld_dexample_1lev \
	$(EXEDIR)/mld_zexample_ml $(EXEDIR)/mld_zexample_1lev \
	$(SMOBJS) $(S1OBJS) $(CMOBJS) $(C1OBJS) \
        $(EXEDIR)/mld_sexample_ml $(EXEDIR)/mld_sexample_1lev \
	$(EXEDIR)/mld_cexample_ml $(EXEDIR)/mld_cexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


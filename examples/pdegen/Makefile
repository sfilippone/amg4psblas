MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.amg4psblas
MLDMODDIR=$(MLDDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDMODDIR) $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).
LINKOPT=

DMOBJS=mld_dexample_ml.o data_input.o mld_dpde_mod.o
D1OBJS=mld_dexample_1lev.o data_input.o mld_dpde_mod.o
SMOBJS=mld_sexample_ml.o data_input.o mld_spde_mod.o
S1OBJS=mld_sexample_1lev.o data_input.o mld_spde_mod.o

EXEDIR=./runs

all: mld_sexample_ml mld_sexample_1lev mld_dexample_ml mld_dexample_1lev


mld_dexample_ml: $(DMOBJS)
	$(FLINK) $(LINKOPT) $(DMOBJS) -o mld_dexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_dexample_ml $(EXEDIR) 

mld_dexample_1lev: $(D1OBJS)
	$(FLINK) $(LINKOPT) $(D1OBJS) -o mld_dexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_dexample_1lev $(EXEDIR) 

mld_dexample_ml.o: data_input.o mld_dpde_mod.o
mld_dexample_1lev.o: data_input.o mld_dpde_mod.o

mld_sexample_ml: $(SMOBJS)
	$(FLINK) $(LINKOPT) $(SMOBJS) -o mld_sexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_sexample_ml $(EXEDIR) 

mld_sexample_1lev: $(S1OBJS)
	$(FLINK) $(LINKOPT) $(S1OBJS) -o mld_sexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_sexample_1lev $(EXEDIR) 

mld_sexample_ml.o: data_input.o  mld_spde_mod.o
mld_sexample_1lev.o: data_input.o mld_spde_mod.o

clean: 
	/bin/rm -f $(DMOBJS) $(D1OBJS) $(SMOBJS) $(S1OBJS) \
        *$(.mod)  $(EXEDIR)/mld_dexample_ml $(EXEDIR)/mld_dexample_1lev\
	 $(EXEDIR)/mld_sexample_ml $(EXEDIR)/mld_sexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


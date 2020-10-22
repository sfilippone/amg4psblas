MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.amg4psblas
MLDMODDIR=$(MLDDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDMODDIR) $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

DFSOBJS=mld_df_sample.o data_input.o
SFSOBJS=mld_sf_sample.o data_input.o
CFSOBJS=mld_cf_sample.o data_input.o
ZFSOBJS=mld_zf_sample.o data_input.o
LINKOPT=

EXEDIR=./runs

all: mld_sf_sample mld_df_sample mld_cf_sample mld_zf_sample 

mld_df_sample: $(DFSOBJS)
	$(FLINK) $(LINKOPT) $(DFSOBJS) -o mld_df_sample \
	$(MLD_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_df_sample $(EXEDIR) 


mld_sf_sample: $(SFSOBJS)
	$(FLINK) $(LINKOPT) $(SFSOBJS) -o mld_sf_sample \
	$(MLD_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_sf_sample $(EXEDIR) 

mld_cf_sample: $(CFSOBJS)
	$(FLINK) $(LINKOPT) $(CFSOBJS) -o mld_cf_sample \
	$(MLD_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_cf_sample $(EXEDIR) 

mld_zf_sample: $(ZFSOBJS)
	$(FLINK) $(LINKOPT) $(ZFSOBJS) -o mld_zf_sample \
	$(MLD_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_zf_sample $(EXEDIR) 

mld_sf_sample.o: data_input.o
mld_df_sample.o: data_input.o
mld_cf_sample.o: data_input.o
mld_zf_sample.o: data_input.o

clean: 
	/bin/rm -f $(DFSOBJS) $(ZFSOBJS) $(SFSOBJS)  $(CFSOBJS) \
        *$(.mod) $(EXEDIR)/mld_sf_sample $(EXEDIR)/mld_cf_sample  \
	$(EXEDIR)/mld_df_sample $(EXEDIR)/mld_zf_sample

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


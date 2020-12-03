AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

DFSOBJS=amg_df_sample.o data_input.o
SFSOBJS=amg_sf_sample.o data_input.o
CFSOBJS=amg_cf_sample.o data_input.o
ZFSOBJS=amg_zf_sample.o data_input.o
LINKOPT=

EXEDIR=./runs

all: amg_sf_sample amg_df_sample amg_cf_sample amg_zf_sample 

amg_df_sample: $(DFSOBJS)
	$(FLINK) $(LINKOPT) $(DFSOBJS) -o amg_df_sample \
	$(AMG_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_df_sample $(EXEDIR) 


amg_sf_sample: $(SFSOBJS)
	$(FLINK) $(LINKOPT) $(SFSOBJS) -o amg_sf_sample \
	$(AMG_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_sf_sample $(EXEDIR) 

amg_cf_sample: $(CFSOBJS)
	$(FLINK) $(LINKOPT) $(CFSOBJS) -o amg_cf_sample \
	$(AMG_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_cf_sample $(EXEDIR) 

amg_zf_sample: $(ZFSOBJS)
	$(FLINK) $(LINKOPT) $(ZFSOBJS) -o amg_zf_sample \
	$(AMG_LIBS)  $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_zf_sample $(EXEDIR) 

amg_sf_sample.o: data_input.o
amg_df_sample.o: data_input.o
amg_cf_sample.o: data_input.o
amg_zf_sample.o: data_input.o

clean: 
	/bin/rm -f $(DFSOBJS) $(ZFSOBJS) $(SFSOBJS)  $(CFSOBJS) \
        *$(.mod) $(EXEDIR)/amg_sf_sample $(EXEDIR)/amg_cf_sample  \
	$(EXEDIR)/amg_df_sample $(EXEDIR)/amg_zf_sample

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)


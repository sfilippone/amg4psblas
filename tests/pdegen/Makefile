AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_krylov -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

LINKOPT=
EXEDIR=./runs

all: amg_s_pde3d amg_d_pde3d amg_s_pde2d amg_d_pde2d 

amg_d_pde3d: amg_d_pde3d.o  data_input.o
	$(FLINK) $(LINKOPT) amg_d_pde3d.o  data_input.o -o amg_d_pde3d $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_d_pde3d $(EXEDIR)

amg_s_pde3d: amg_s_pde3d.o   data_input.o
	$(FLINK) $(LINKOPT) amg_s_pde3d.o  data_input.o -o amg_s_pde3d $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_s_pde3d $(EXEDIR)

amg_d_pde2d: amg_d_pde2d.o  data_input.o
	$(FLINK) $(LINKOPT) amg_d_pde2d.o  data_input.o -o amg_d_pde2d $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_d_pde2d $(EXEDIR)

amg_s_pde2d: amg_s_pde2d.o   data_input.o
	$(FLINK) $(LINKOPT) amg_s_pde2d.o  data_input.o -o amg_s_pde2d $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_s_pde2d $(EXEDIR)

amg_d_pde3d_rebld: amg_d_pde3d_rebld.o  data_input.o
	$(FLINK) $(LINKOPT) amg_d_pde3d_rebld.o  data_input.o -o amg_d_pde3d_rebld $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_d_pde3d_rebld $(EXEDIR)

amg_d_pde3d.o amg_s_pde3d.o amg_d_pde2d.o amg_s_pde2d.o: data_input.o 

check: all
	cd runs && ./amg_d_pde2d <amg_pde2d.inp && ./amg_s_pde2d<amg_pde2d.inp


clean: 
	/bin/rm -f data_input.o amg_d_pde3d.o amg_s_pde3d.o amg_d_pde2d.o amg_s_pde2d.o *$(.mod)\
        $(EXEDIR)/mld_d_pde3d  $(EXEDIR)/mld_s_pde3d $(EXEDIR)/mld_d_pde2d  $(EXEDIR)/mld_s_pde2d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)




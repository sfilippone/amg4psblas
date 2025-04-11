AMGDIR=../..
AMGINCDIR=$(AMGDIR)/include
include $(AMGINCDIR)/Make.inc.amg4psblas
AMGMODDIR=$(AMGDIR)/modules
AMGLIBDIR=$(AMGDIR)/lib
AMG_LIBS=-L$(AMGLIBDIR) -lpsb_linsolve -lamg_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(AMGMODDIR) $(FMFLAG)$(AMGINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

PD3DOBJS=amg_pde3d_newslv.o  data_input.o amg_d_tlu_solver.o amg_d_tlu_solver_impl.o
PSOBJS=spde.o  data_input.o 
EXEDIR=./runs
LINKOPT=

all: amg_pde3d_newslv 


amg_pde3d_newslv: $(PD3DOBJS)
	$(FLINK) $(LINKOPT) $(PD3DOBJS) -o amg_pde3d_newslv $(AMG_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv amg_pde3d_newslv $(EXEDIR)


amg_pde3d_newslv.o amg_d_tlu_solver_impl.o: data_input.o amg_d_tlu_solver.o


clean: 
	/bin/rm -f $(PD3DOBJS) $(EXEDIR)/amg_pde3d_newslv 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)




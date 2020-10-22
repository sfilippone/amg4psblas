MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.amg4psblas
MLDMODDIR=$(MLDDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDMODDIR) $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

PD3DOBJS=mld_pde3d_newslv.o  data_input.o mld_d_tlu_solver.o mld_d_tlu_solver_impl.o
PSOBJS=spde.o  data_input.o 
EXEDIR=./runs
LINKOPT=

all: mld_pde3d_newslv 


mld_pde3d_newslv: $(PD3DOBJS)
	$(FLINK) $(LINKOPT) $(PD3DOBJS) -o mld_pde3d_newslv $(MLD_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_pde3d_newslv $(EXEDIR)


mld_pde3d_newslv.o mld_d_tlu_solver_impl.o: data_input.o mld_d_tlu_solver.o


clean: 
	/bin/rm -f $(PD3DOBJS) $(EXEDIR)/mld_pde3d_newslv 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)




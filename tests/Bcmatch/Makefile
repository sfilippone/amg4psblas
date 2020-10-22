MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.amg4psblas
MLDMODDIR=$(MLDDIR)/modules
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDMODDIR) $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

HSL_DIR=/opt/hsl/2.3.1/gnu/6.4.0
HSL_INCDIR=$(HSL_DIR)/include
HSL_LIBDIR=$(HSL_DIR)/lib
HSL_LIBS=-lhsl_mc64 -L$(HSL_LIBDIR)
HSL_FLAGS= -DHAVE_HSL -I$(HSL_INCDIR)

# SPRAL package for auction algorithm
SPRAL_DIR=/opt/spral/2015.04.20/gnu/6.4.0
SPRAL_INCDIR=$(SPRAL_DIR)/include
SPRAL_LIBDIR=$(SPRAL_DIR)/lib
SPRAL_LIBS=-lspral -L$(SPRAL_LIBDIR)
SPRAL_FLAGS=-DHAVE_SPRAL -I$(SPRAL_INCDIR)

BCM_DIR=/opt/bcm/0.9/gnu/6.4.0
BCM_INCDIR=$(BCM_DIR)/include
BCM_LIBDIR=$(BCM_DIR)/lib
BCM_LDLIBS=-lBCM -L$(BCM_LIBDIR) $(HSL_LIBS) $(SPRAL_LIBS)

CDEFINES=$(MLDCDEFINES) -I$(BCM_INCDIR)

LINKOPT=
EXEDIR=./runs

all:  mld_d_pde3d 

BCMOBJS= mld_d_bcmatch_aggregator_mod.o mld_d_bcmatch_aggregator_mat_asb.o \
         mld_d_bcmatch_aggregator_tprol.o  mld_d_bcmatch_map_to_tprol.o \
         mld_daggrmat_unsmth_spmm_asb.o bootCMatch_interface.o

mld_d_pde3d: mld_d_pde3d.o  data_input.o $(BCMOBJS)
	$(FLINK) $(LINKOPT) mld_d_pde3d.o  data_input.o $(BCMOBJS) -o mld_d_pde3d $(MLD_LIBS) $(BCM_LDLIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv mld_d_pde3d $(EXEDIR)

mld_d_pde3d.o: data_input.o mld_d_bcmatch_aggregator_mod.o
mld_d_bcmatch_aggregator_mat_asb.o mld_d_bcmatch_aggregator_tprol.o  mld_d_bcmatch_map_to_tprol.o: mld_d_bcmatch_aggregator_mod.o

check: all
	cd runs && ./mld_d_pde2d <mld_pde2d.inp && ./mld_s_pde2d<mld_pde2d.inp


clean: 
	/bin/rm -f data_input.o mld_d_pde3d.o *$(.mod)  $(BCMOBJS)\
        $(EXEDIR)/mld_d_pde3d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)




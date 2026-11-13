AMG_DIR ?= $(abspath ../../..)
ifneq ($(wildcard $(AMG_DIR)/Make.inc),)
AMG_MAKEINC ?= $(AMG_DIR)/Make.inc
include $(AMG_MAKEINC)
AMG_TEST_INCDIR=$(AMG_DIR)/include
AMG_TEST_LIBDIR=$(AMG_DIR)/lib
else
AMG_MAKEINC ?= $(AMG_DIR)/include/Make.inc.amg4psblas
include $(AMG_MAKEINC)
AMG_TEST_INCDIR=$(INSTALL_INCLUDEDIR)
AMG_TEST_LIBDIR=$(INSTALL_LIBDIR)
endif

CINCLUDES=-I$(AMG_TEST_INCDIR) -I$(PSBLAS_INCDIR)
AMGC_LIBS=-L$(AMG_TEST_LIBDIR) -lamg_cbind -lamg_prec
PSBC_LIBS=-L$(PSBLAS_LIBDIR) -lpsb_cbind -lpsb_linsolve -lpsb_prec -lpsb_ext
HEADERS=$(wildcard $(AMG_TEST_INCDIR)/*.h $(PSBLAS_INCDIR)/*.h)
LIBRARIES=$(AMG_TEST_LIBDIR)/libamg_cbind.a $(AMG_TEST_LIBDIR)/libamg_prec.a \
	$(PSBLAS_LIBDIR)/libpsb_cbind.a $(PSBLAS_LIBDIR)/libpsb_linsolve.a \
	$(PSBLAS_LIBDIR)/libpsb_prec.a $(PSBLAS_LIBDIR)/libpsb_ext.a \
	$(PSBLAS_LIBDIR)/libpsb_util.a $(PSBLAS_LIBDIR)/libpsb_base.a
MPIEXEC ?= mpiexec
MPIEXEC_NUMPROC_FLAG ?= -n
MPIEXEC_PREFLAGS ?=
MPIEXEC_POSTFLAGS ?=
NP ?= 2
LINKOPT ?=

.PHONY: all check clean
all: prec_extensions

prec_extensions.o: prec_extensions.c Makefile $(AMG_MAKEINC) \
	$(wildcard $(PSBLAS_INCDIR)/Make.inc.psblas) $(HEADERS)
	$(MPCC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -std=c99 -c $<

prec_extensions: prec_extensions.o $(LIBRARIES)
	$(FLINK) $(LINKOPT) -o $@ $< $(AMGC_LIBS) $(PSBC_LIBS) \
	  $(PSBCLDLIBS) $(PSBLAS_LIBS) $(PSBGPULDLIBS) $(LDLIBS) -lm

check: prec_extensions
	$(MPIEXEC) $(MPIEXEC_NUMPROC_FLAG) $(NP) $(MPIEXEC_PREFLAGS) \
	  ./prec_extensions $(MPIEXEC_POSTFLAGS)

clean:
	rm -f prec_extensions prec_extensions.o

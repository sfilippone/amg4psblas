TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=.
DEST=../

CINCLUDES=-I. -I$(INCDIR) -I$(PSBLAS_INCDIR)
FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(PSBLAS_INCLUDES)
AMGFDEFINES=-DPSB_HAVE_LAPACK -DPSB_HAVE_FLUSH_STMT -DPSB_MPI_MOD  $(PSBFDEFINES) $(FCUDEFINES)
FDEFINES=$(AMGFDEFINES) 


OBJS=amg_prec_cbind_mod.o amg_dprec_cbind_mod.o  amg_c_dprec.o amg_zprec_cbind_mod.o  amg_c_zprec.o amg_cbind_ext_mod.o
CMOD=amg_cbind.h amg_c_dprec.h amg_c_zprec.h amg_c_prec_info.h


LIBMOD=amg_prec_cbind_mod$(.mod) amg_dprec_cbind_mod$(.mod) amg_zprec_cbind_mod$(.mod) amg_cbind_ext_mod$(.mod)
LOCAL_MODS=$(LIBMOD)
#LIBNAME=$(CPRECLIBNAME)

.PHONY: objs lib clean veryclean

objs: $(OBJS) $(CMOD)
	/bin/cp -p $(LIBMOD) $(MODDIR)
	/bin/cp -p $(CMOD) $(INCDIR)

lib: objs $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(DEST)

amg_prec_cbind_mod.o: amg_dprec_cbind_mod.o amg_zprec_cbind_mod.o amg_cbind_ext_mod.o

$(OBJS): Makefile $(TOP)/Make.inc $(wildcard $(PSBLAS_INCDIR)/Make.inc.psblas)

amg_c_dprec.o amg_c_zprec.o: $(CMOD) $(INCDIR)/amg_config.h \
	$(wildcard $(PSBLAS_INCDIR)/*.h)

amg_dprec_cbind_mod.o amg_zprec_cbind_mod.o amg_cbind_ext_mod.o: \
	$(MODDIR)/amg_prec_mod$(.mod)

amg_prec_cbind_mod.o amg_dprec_cbind_mod.o amg_zprec_cbind_mod.o amg_cbind_ext_mod.o: \
	$(PSBLAS_MODDIR)/psb_base_cbind_mod$(.mod)

amg_prec_cbind_mod.o amg_dprec_cbind_mod.o amg_zprec_cbind_mod.o: \
	$(PSBLAS_MODDIR)/psb_prec_cbind_mod$(.mod)

amg_dprec_cbind_mod.o amg_zprec_cbind_mod.o: \
	$(PSBLAS_MODDIR)/psb_base_mod$(.mod) \
	$(PSBLAS_MODDIR)/psb_prec_mod$(.mod) \
	$(PSBLAS_MODDIR)/psb_linsolve_mod$(.mod) \
	$(PSBLAS_MODDIR)/psb_objhandle_mod$(.mod) \
	$(wildcard $(PSBLAS_MODDIR)/psb_cuda_mod$(.mod))

amg_dprec_cbind_mod.o: $(PSBLAS_MODDIR)/psb_dlinsolve_cbind_mod$(.mod)
amg_zprec_cbind_mod.o: $(PSBLAS_MODDIR)/psb_zlinsolve_cbind_mod$(.mod)

amg_cbind_ext_mod.o: $(MODDIR)/amg_d_jac_smoother$(.mod) \
	$(MODDIR)/amg_z_jac_smoother$(.mod)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean

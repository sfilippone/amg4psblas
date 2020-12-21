TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=.
DEST=../

CINCLUDES=-I. -I$(LIBDIR) -I$(PSBLAS_INCDIR)
FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(PSBLAS_INCLUDES)


OBJS=amg_prec_cbind_mod.o amg_dprec_cbind_mod.o  amg_c_dprec.o amg_zprec_cbind_mod.o  amg_c_zprec.o
CMOD=amg_cbind.h amg_c_dprec.h  amg_c_zprec.h amg_const.h


LIBMOD=amg_prec_cbind_mod$(.mod) amg_dprec_cbind_mod$(.mod) amg_zprec_cbind_mod$(.mod)
LOCAL_MODS=$(LIBMOD)
#LIBNAME=$(CPRECLIBNAME)


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(DEST)
	/bin/cp -p $(LIBMOD) $(CMOD) $(DEST)

amg_prec_cbind_mod.o: amg_dprec_cbind_mod.o  amg_zprec_cbind_mod.o
#amg_prec_cbind_mod.o: psb_prec_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean

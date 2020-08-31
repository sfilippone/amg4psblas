include ../../Make.inc
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules 
HERE=..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)



DMPFOBJS=

SMPFOBJS=

ZMPFOBJS=

CMPFOBJS=


MPFOBJS=$(SMPFOBJS) $(DMPFOBJS) $(CMPFOBJS) $(ZMPFOBJS) 

MPCOBJS=amg_dslud_interface.o amg_zslud_interface.o


DINNEROBJS= amg_dmlprec_bld.o amg_dfile_prec_descr.o \
	amg_d_smoothers_bld.o amg_d_hierarchy_bld.o \
	amg_dmlprec_aply.o   \
	$(DMPFOBJS) amg_d_extprol_bld.o

SINNEROBJS=   amg_smlprec_bld.o amg_sfile_prec_descr.o \
	amg_s_smoothers_bld.o amg_s_hierarchy_bld.o \
	amg_smlprec_aply.o   \
	$(SMPFOBJS) amg_s_extprol_bld.o 

ZINNEROBJS=   amg_zmlprec_bld.o amg_zfile_prec_descr.o \
	amg_z_smoothers_bld.o amg_z_hierarchy_bld.o \
	amg_zmlprec_aply.o   \
	$(ZMPFOBJS)  amg_z_extprol_bld.o

CINNEROBJS=   amg_cmlprec_bld.o amg_cfile_prec_descr.o \
	amg_c_smoothers_bld.o amg_c_hierarchy_bld.o \
	amg_cmlprec_aply.o   \
	$(CMPFOBJS)  amg_c_extprol_bld.o

INNEROBJS= $(SINNEROBJS) $(DINNEROBJS) $(CINNEROBJS) $(ZINNEROBJS)


DOUTEROBJS=amg_dprecbld.o  amg_dprecset.o amg_dprecinit.o amg_dprecaply.o amg_dcprecset.o

SOUTEROBJS=amg_sprecbld.o  amg_sprecset.o amg_sprecinit.o amg_sprecaply.o amg_scprecset.o

ZOUTEROBJS=amg_zprecbld.o  amg_zprecset.o amg_zprecinit.o amg_zprecaply.o amg_zcprecset.o

COUTEROBJS=amg_cprecbld.o  amg_cprecset.o amg_cprecinit.o amg_cprecaply.o amg_ccprecset.o


OUTEROBJS=$(SOUTEROBJS) $(DOUTEROBJS) $(COUTEROBJS) $(ZOUTEROBJS)

F90OBJS=$(OUTEROBJS) $(INNEROBJS)

COBJS= amg_sslu_interface.o  \
       amg_dslu_interface.o amg_dumf_interface.o \
       amg_cslu_interface.o  \
       amg_zslu_interface.o amg_zumf_interface.o

OBJS=$(F90OBJS) $(COBJS) $(MPCOBJS)

LIBNAME=libamg_prec.a

lib:  $(OBJS) aggrd levd smoothd solvd
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

aggrd:
	$(MAKE) -C aggregator
levd:
	$(MAKE) -C level
smoothd:
	$(MAKE) -C smoother
solvd:
	$(MAKE) -C solver

mpobjs: 
	(make $(MPFOBJS) FC="$(MPFC)" FCOPT="$(FCOPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: solvclean smoothclean levclean aggrclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

aggrclean:
		$(MAKE) -C aggregator clean
levclean:
		$(MAKE) -C level clean
smoothclean:
		$(MAKE) -C smoother clean
solvclean:
		$(MAKE) -C solver clean

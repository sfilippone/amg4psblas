include ../../../Make.inc
LIBDIR=../../../lib
INCDIR=../../../include
MODDIR=../../../modules 
HERE=../..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)



OBJS=amg_c_base_onelev_build.o \
amg_c_base_onelev_check.o \
amg_c_base_onelev_cnv.o \
amg_c_base_onelev_csetc.o \
amg_c_base_onelev_cseti.o \
amg_c_base_onelev_csetr.o \
amg_c_base_onelev_descr.o \
amg_c_base_onelev_dump.o \
amg_c_base_onelev_free.o \
amg_c_base_onelev_mat_asb.o \
amg_c_base_onelev_setag.o \
amg_c_base_onelev_setsm.o \
amg_c_base_onelev_setsv.o \
amg_c_base_onelev_map_rstr.o \
amg_c_base_onelev_map_prol.o \
amg_d_base_onelev_build.o \
amg_d_base_onelev_check.o \
amg_d_base_onelev_cnv.o \
amg_d_base_onelev_csetc.o \
amg_d_base_onelev_cseti.o \
amg_d_base_onelev_csetr.o \
amg_d_base_onelev_descr.o \
amg_d_base_onelev_dump.o \
amg_d_base_onelev_free.o \
amg_d_base_onelev_mat_asb.o \
amg_d_base_onelev_setag.o \
amg_d_base_onelev_setsm.o \
amg_d_base_onelev_setsv.o \
amg_d_base_onelev_map_rstr.o \
amg_d_base_onelev_map_prol.o \
amg_s_base_onelev_build.o \
amg_s_base_onelev_check.o \
amg_s_base_onelev_cnv.o \
amg_s_base_onelev_csetc.o \
amg_s_base_onelev_cseti.o \
amg_s_base_onelev_csetr.o \
amg_s_base_onelev_descr.o \
amg_s_base_onelev_dump.o \
amg_s_base_onelev_free.o \
amg_s_base_onelev_mat_asb.o \
amg_s_base_onelev_setag.o \
amg_s_base_onelev_setsm.o \
amg_s_base_onelev_setsv.o \
amg_s_base_onelev_map_rstr.o \
amg_s_base_onelev_map_prol.o \
amg_z_base_onelev_build.o \
amg_z_base_onelev_check.o \
amg_z_base_onelev_cnv.o \
amg_z_base_onelev_csetc.o \
amg_z_base_onelev_cseti.o \
amg_z_base_onelev_csetr.o \
amg_z_base_onelev_descr.o \
amg_z_base_onelev_dump.o \
amg_z_base_onelev_free.o \
amg_z_base_onelev_mat_asb.o \
amg_z_base_onelev_setag.o \
amg_z_base_onelev_setsm.o \
amg_z_base_onelev_setsv.o \
amg_z_base_onelev_map_rstr.o \
amg_z_base_onelev_map_prol.o 


LIBNAME=libamg_prec.a

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

mpobjs: 
	(make $(MPFOBJS) FC="$(MPFC)" FCOPT="$(FCOPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

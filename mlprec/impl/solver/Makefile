include ../../../Make.inc
LIBDIR=../../../lib
INCDIR=../../../include
MODDIR=../../../modules 
HERE=../..
FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)



OBJS=mld_c_base_solver_apply.o \
mld_c_base_solver_apply_vect.o \
mld_c_base_solver_bld.o \
mld_c_base_solver_check.o \
mld_c_base_solver_clone.o \
mld_c_base_solver_clone_settings.o \
mld_c_base_solver_clear_data.o \
mld_c_base_solver_cnv.o \
mld_c_base_solver_csetc.o \
mld_c_base_solver_cseti.o \
mld_c_base_solver_csetr.o \
mld_c_base_solver_descr.o \
mld_c_base_solver_dmp.o \
mld_c_base_solver_free.o \
mld_c_diag_solver_apply.o \
mld_c_diag_solver_apply_vect.o \
mld_c_diag_solver_bld.o \
mld_c_diag_solver_dmp.o \
mld_c_diag_solver_clone.o \
mld_c_diag_solver_clear_data.o \
mld_c_diag_solver_cnv.o \
mld_c_gs_solver_bld.o \
mld_c_gs_solver_clone.o \
mld_c_gs_solver_cnv.o \
mld_c_gs_solver_dmp.o \
mld_c_gs_solver_apply.o \
mld_c_gs_solver_apply_vect.o \
mld_c_gs_solver_clear_data.o \
mld_c_gs_solver_clone_settings.o \
mld_c_bwgs_solver_bld.o \
mld_c_bwgs_solver_apply.o \
mld_c_bwgs_solver_apply_vect.o \
mld_c_id_solver_apply.o \
mld_c_id_solver_apply_vect.o \
mld_c_id_solver_clone.o \
mld_c_ilu_solver_apply.o \
mld_c_ilu_solver_apply_vect.o \
mld_c_ilu_solver_bld.o \
mld_c_ilu_solver_clone.o \
mld_c_ilu_solver_clear_data.o \
mld_c_ilu_solver_clone_settings.o \
mld_c_ilu_solver_cnv.o \
mld_c_ilu_solver_dmp.o \
mld_c_mumps_solver_apply.o \
mld_c_mumps_solver_apply_vect.o \
mld_c_mumps_solver_bld.o \
mld_d_base_solver_apply.o \
mld_d_base_solver_apply_vect.o \
mld_d_base_solver_bld.o \
mld_d_base_solver_check.o \
mld_d_base_solver_clone.o \
mld_d_base_solver_clone_settings.o \
mld_d_base_solver_clear_data.o \
mld_d_base_solver_cnv.o \
mld_d_base_solver_csetc.o \
mld_d_base_solver_cseti.o \
mld_d_base_solver_csetr.o \
mld_d_base_solver_descr.o \
mld_d_base_solver_dmp.o \
mld_d_base_solver_free.o \
mld_d_diag_solver_apply.o \
mld_d_diag_solver_apply_vect.o \
mld_d_diag_solver_bld.o \
mld_d_diag_solver_dmp.o \
mld_d_diag_solver_clone.o \
mld_d_diag_solver_clear_data.o \
mld_d_diag_solver_cnv.o \
mld_d_bwgs_solver_bld.o \
mld_d_bwgs_solver_apply.o \
mld_d_bwgs_solver_apply_vect.o \
mld_d_gs_solver_bld.o \
mld_d_gs_solver_clone.o \
mld_d_gs_solver_cnv.o \
mld_d_gs_solver_dmp.o \
mld_d_gs_solver_apply.o \
mld_d_gs_solver_apply_vect.o \
mld_d_gs_solver_clear_data.o \
mld_d_gs_solver_clone_settings.o \
mld_d_id_solver_apply.o \
mld_d_id_solver_apply_vect.o \
mld_d_id_solver_clone.o \
mld_d_ilu_solver_apply.o \
mld_d_ilu_solver_apply_vect.o \
mld_d_ilu_solver_bld.o \
mld_d_ilu_solver_clone.o \
mld_d_ilu_solver_clear_data.o \
mld_d_ilu_solver_clone_settings.o \
mld_d_ilu_solver_cnv.o \
mld_d_ilu_solver_dmp.o \
mld_d_mumps_solver_apply.o \
mld_d_mumps_solver_apply_vect.o \
mld_d_mumps_solver_bld.o \
mld_s_base_solver_apply.o \
mld_s_base_solver_apply_vect.o \
mld_s_base_solver_bld.o \
mld_s_base_solver_check.o \
mld_s_base_solver_clone.o \
mld_s_base_solver_clone_settings.o \
mld_s_base_solver_clear_data.o \
mld_s_base_solver_cnv.o \
mld_s_base_solver_csetc.o \
mld_s_base_solver_cseti.o \
mld_s_base_solver_csetr.o \
mld_s_base_solver_descr.o \
mld_s_base_solver_dmp.o \
mld_s_base_solver_free.o \
mld_s_diag_solver_apply.o \
mld_s_diag_solver_apply_vect.o \
mld_s_diag_solver_bld.o \
mld_s_diag_solver_dmp.o \
mld_s_diag_solver_clone.o \
mld_s_diag_solver_clear_data.o \
mld_s_diag_solver_cnv.o \
mld_s_gs_solver_bld.o \
mld_s_gs_solver_clone.o \
mld_s_gs_solver_cnv.o \
mld_s_gs_solver_dmp.o \
mld_s_gs_solver_apply.o \
mld_s_gs_solver_apply_vect.o \
mld_s_gs_solver_clear_data.o \
mld_s_gs_solver_clone_settings.o \
mld_s_bwgs_solver_bld.o \
mld_s_bwgs_solver_apply.o \
mld_s_bwgs_solver_apply_vect.o \
mld_s_id_solver_apply.o \
mld_s_id_solver_apply_vect.o \
mld_s_id_solver_clone.o \
mld_s_ilu_solver_apply.o \
mld_s_ilu_solver_apply_vect.o \
mld_s_ilu_solver_bld.o \
mld_s_ilu_solver_clone.o \
mld_s_ilu_solver_clear_data.o \
mld_s_ilu_solver_clone_settings.o \
mld_s_ilu_solver_cnv.o \
mld_s_ilu_solver_dmp.o \
mld_s_mumps_solver_apply.o \
mld_s_mumps_solver_apply_vect.o \
mld_s_mumps_solver_bld.o \
mld_z_base_solver_apply.o \
mld_z_base_solver_apply_vect.o \
mld_z_base_solver_bld.o \
mld_z_base_solver_check.o \
mld_z_base_solver_clone.o \
mld_z_base_solver_clone_settings.o \
mld_z_base_solver_clear_data.o \
mld_z_base_solver_cnv.o \
mld_z_base_solver_csetc.o \
mld_z_base_solver_cseti.o \
mld_z_base_solver_csetr.o \
mld_z_base_solver_descr.o \
mld_z_base_solver_dmp.o \
mld_z_base_solver_free.o \
mld_z_diag_solver_apply.o \
mld_z_diag_solver_apply_vect.o \
mld_z_diag_solver_bld.o \
mld_z_diag_solver_dmp.o \
mld_z_diag_solver_clone.o \
mld_z_diag_solver_clear_data.o \
mld_z_diag_solver_cnv.o \
mld_z_gs_solver_bld.o \
mld_z_gs_solver_clone.o \
mld_z_gs_solver_cnv.o \
mld_z_gs_solver_dmp.o \
mld_z_gs_solver_apply.o \
mld_z_gs_solver_apply_vect.o \
mld_z_gs_solver_clear_data.o \
mld_z_gs_solver_clone_settings.o \
mld_z_bwgs_solver_bld.o \
mld_z_bwgs_solver_apply.o \
mld_z_bwgs_solver_apply_vect.o \
mld_z_id_solver_apply.o \
mld_z_id_solver_apply_vect.o \
mld_z_id_solver_clone.o \
mld_z_ilu_solver_clear_data.o \
mld_z_ilu_solver_clone_settings.o \
mld_z_ilu_solver_apply.o \
mld_z_ilu_solver_apply_vect.o \
mld_z_ilu_solver_bld.o \
mld_z_ilu_solver_clone.o \
mld_z_ilu_solver_cnv.o \
mld_z_ilu_solver_dmp.o \
mld_z_mumps_solver_apply.o \
mld_z_mumps_solver_apply_vect.o \
mld_z_mumps_solver_bld.o \


LIBNAME=libmld_prec.a

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


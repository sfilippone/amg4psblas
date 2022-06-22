include ../Make.inc
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
HERE=.

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)


DMODOBJS=amg_d_prec_type.o \
	amg_d_inner_mod.o amg_d_ilu_solver.o amg_d_diag_solver.o amg_d_jac_smoother.o amg_d_as_smoother.o \
	amg_d_umf_solver.o amg_d_slu_solver.o amg_d_sludist_solver.o amg_d_id_solver.o\
	amg_d_base_solver_mod.o amg_d_base_smoother_mod.o amg_d_onelev_mod.o \
	amg_d_gs_solver.o amg_d_mumps_solver.o \
	amg_d_base_aggregator_mod.o \
	amg_d_dec_aggregator_mod.o amg_d_symdec_aggregator_mod.o \
	amg_d_ainv_solver.o amg_d_base_ainv_mod.o \
	amg_d_invk_solver.o amg_d_invt_solver.o amg_d_krm_solver.o \
	amg_d_matchboxp_mod.o amg_d_parmatch_aggregator_mod.o

SMODOBJS=amg_s_prec_type.o amg_s_ilu_fact_mod.o \
	amg_s_inner_mod.o amg_s_ilu_solver.o amg_s_diag_solver.o amg_s_jac_smoother.o amg_s_as_smoother.o \
	amg_s_slu_solver.o amg_s_id_solver.o\
	amg_s_base_solver_mod.o amg_s_base_smoother_mod.o amg_s_onelev_mod.o  \
	amg_s_gs_solver.o amg_s_mumps_solver.o \
	amg_s_base_aggregator_mod.o  \
	amg_s_dec_aggregator_mod.o amg_s_symdec_aggregator_mod.o \
	amg_s_ainv_solver.o amg_s_base_ainv_mod.o \
	amg_s_invk_solver.o amg_s_invt_solver.o amg_s_krm_solver.o \
	amg_s_matchboxp_mod.o amg_s_parmatch_aggregator_mod.o

ZMODOBJS=amg_z_prec_type.o amg_z_ilu_fact_mod.o \
	amg_z_inner_mod.o amg_z_ilu_solver.o amg_z_diag_solver.o amg_z_jac_smoother.o amg_z_as_smoother.o \
	amg_z_umf_solver.o amg_z_slu_solver.o amg_z_sludist_solver.o amg_z_id_solver.o\
	amg_z_base_solver_mod.o amg_z_base_smoother_mod.o amg_z_onelev_mod.o \
	amg_z_gs_solver.o amg_z_mumps_solver.o \
	amg_z_base_aggregator_mod.o  \
	amg_z_dec_aggregator_mod.o amg_z_symdec_aggregator_mod.o \
	amg_z_ainv_solver.o amg_z_base_ainv_mod.o \
	amg_z_invk_solver.o amg_z_invt_solver.o amg_z_krm_solver.o

CMODOBJS=amg_c_prec_type.o amg_c_ilu_fact_mod.o \
	amg_c_inner_mod.o amg_c_ilu_solver.o amg_c_diag_solver.o amg_c_jac_smoother.o amg_c_as_smoother.o \
	amg_c_slu_solver.o amg_c_id_solver.o\
	amg_c_base_solver_mod.o amg_c_base_smoother_mod.o amg_c_onelev_mod.o  \
	amg_c_gs_solver.o amg_c_mumps_solver.o \
	amg_c_base_aggregator_mod.o  \
	amg_c_dec_aggregator_mod.o amg_c_symdec_aggregator_mod.o \
	amg_c_ainv_solver.o amg_c_base_ainv_mod.o \
	amg_c_invk_solver.o amg_c_invt_solver.o amg_c_krm_solver.o



MODOBJS=amg_base_prec_type.o amg_prec_type.o  amg_prec_mod.o \
	 amg_s_prec_mod.o  amg_d_prec_mod.o amg_c_prec_mod.o amg_z_prec_mod.o \
	 amg_base_ainv_mod.o amg_ainv_mod.o\
	$(SMODOBJS) $(DMODOBJS) $(CMODOBJS) $(ZMODOBJS)


OBJS=$(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=libamg_prec.a

all: objs impld

objs: $(OBJS)
	/bin/cp -p amg_const.h $(INCDIR)
	/bin/cp -p *$(.mod) $(MODDIR)

impld: objs
	cd impl && $(MAKE) 

lib: $(OBJS) impld
	cd impl && $(MAKE) lib
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)


$(MODOBJS): $(PSBLAS_MODDIR)/$(PSBBASEMODNAME)$(.mod)

amg_base_prec_type.o: amg_const.h
amg_s_prec_type.o amg_d_prec_type.o amg_c_prec_type.o amg_z_prec_type.o	: amg_base_prec_type.o
amg_prec_type.o: amg_s_prec_type.o amg_d_prec_type.o amg_c_prec_type.o amg_z_prec_type.o
amg_prec_mod.o: amg_prec_type.o amg_s_prec_mod.o amg_d_prec_mod.o amg_c_prec_mod.o amg_z_prec_mod.o
amg_s_krm_solver.o: amg_s_prec_type.o  amg_s_base_solver_mod.o
amg_d_krm_solver.o: amg_d_prec_type.o  amg_d_base_solver_mod.o
amg_c_krm_solver.o: amg_c_prec_type.o  amg_c_base_solver_mod.o
amg_z_krm_solver.o: amg_z_prec_type.o  amg_z_base_solver_mod.o

amg_s_prec_mod.o: amg_s_krm_solver.o
amg_d_prec_mod.o: amg_d_krm_solver.o
amg_c_prec_mod.o: amg_c_krm_solver.o
amg_z_prec_mod.o: amg_z_krm_solver.o

$(SINNEROBJS) $(SOUTEROBJS): $(SMODOBJS)
$(DINNEROBJS) $(DOUTEROBJS): $(DMODOBJS)
$(CINNEROBJS) $(COUTEROBJS): $(CMODOBJS)
$(ZINNEROBJS) $(ZOUTEROBJS): $(ZMODOBJS)

amg_s_inner_mod.o: amg_s_prec_type.o
amg_d_inner_mod.o: amg_d_prec_type.o
amg_c_inner_mod.o: amg_c_prec_type.o
amg_z_inner_mod.o: amg_z_prec_type.o

amg_s_prec_mod.o: $(SMODOBJS)
amg_d_prec_mod.o: $(DMODOBJS)
amg_c_prec_mod.o: $(CMODOBJS)
amg_z_prec_mod.o: $(ZMODOBJS)


amg_s_prec_type.o: amg_s_onelev_mod.o
amg_d_prec_type.o: amg_d_onelev_mod.o
amg_c_prec_type.o: amg_c_onelev_mod.o
amg_z_prec_type.o: amg_z_onelev_mod.o

amg_s_onelev_mod.o: amg_s_base_smoother_mod.o amg_s_dec_aggregator_mod.o amg_s_parmatch_aggregator_mod.o 
amg_d_onelev_mod.o: amg_d_base_smoother_mod.o amg_d_dec_aggregator_mod.o amg_d_parmatch_aggregator_mod.o 
amg_c_onelev_mod.o: amg_c_base_smoother_mod.o amg_c_dec_aggregator_mod.o
amg_z_onelev_mod.o: amg_z_base_smoother_mod.o amg_z_dec_aggregator_mod.o

amg_s_base_aggregator_mod.o: amg_base_prec_type.o
amg_s_parmatch_aggregator_mod.o amg_s_dec_aggregator_mod.o: amg_s_base_aggregator_mod.o
amg_s_hybrid_aggregator_mod.o amg_s_symdec_aggregator_mod.o: amg_s_dec_aggregator_mod.o
amg_s_parmatch_aggregator_mod.o: amg_s_matchboxp_mod.o 

amg_d_base_aggregator_mod.o: amg_base_prec_type.o
amg_d_parmatch_aggregator_mod.o amg_d_dec_aggregator_mod.o: amg_d_base_aggregator_mod.o
amg_d_hybrid_aggregator_mod.o amg_d_symdec_aggregator_mod.o: amg_d_dec_aggregator_mod.o
amg_d_parmatch_aggregator_mod.o: amg_d_matchboxp_mod.o 

amg_c_base_aggregator_mod.o: amg_base_prec_type.o
amg_c_parmatch_aggregator_mod.o amg_c_dec_aggregator_mod.o: amg_c_base_aggregator_mod.o
amg_c_hybrid_aggregator_mod.o amg_c_symdec_aggregator_mod.o: amg_c_dec_aggregator_mod.o

amg_z_base_aggregator_mod.o: amg_base_prec_type.o
amg_z_parmatch_aggregator_mod.o amg_z_dec_aggregator_mod.o: amg_z_base_aggregator_mod.o
amg_z_hybrid_aggregator_mod.o amg_z_symdec_aggregator_mod.o: amg_z_dec_aggregator_mod.o

amg_s_base_smoother_mod.o: amg_s_base_solver_mod.o
amg_d_base_smoother_mod.o: amg_d_base_solver_mod.o
amg_c_base_smoother_mod.o: amg_c_base_solver_mod.o
amg_z_base_smoother_mod.o: amg_z_base_solver_mod.o

amg_s_ainv_solver.o: amg_s_base_ainv_mod.o
amg_c_ainv_solver.o: amg_c_base_ainv_mod.o
amg_d_ainv_solver.o: amg_d_base_ainv_mod.o
amg_z_ainv_solver.o: amg_z_base_ainv_mod.o
amg_s_base_ainv_mod.o: amg_s_base_solver_mod.o amg_base_ainv_mod.o
amg_c_base_ainv_mod.o: amg_c_base_solver_mod.o amg_base_ainv_mod.o
amg_d_base_ainv_mod.o: amg_d_base_solver_mod.o amg_base_ainv_mod.o
amg_z_base_ainv_mod.o: amg_z_base_solver_mod.o amg_base_ainv_mod.o

amg_s_base_solver_mod.o amg_d_base_solver_mod.o amg_c_base_solver_mod.o amg_z_base_solver_mod.o: amg_base_prec_type.o

amg_d_mumps_solver.o amg_d_gs_solver.o amg_d_id_solver.o amg_d_sludist_solver.o amg_d_slu_solver.o \
amg_d_umf_solver.o amg_d_diag_solver.o amg_d_ilu_solver.o: amg_d_base_solver_mod.o  amg_d_prec_type.o

#amg_d_ilu_fact_mod.o: amg_base_prec_type.o  amg_d_base_solver_mod.o
#amg_d_ilu_solver.o amg_d_iluk_fact.o: amg_d_ilu_fact_mod.o
amg_d_as_smoother.o amg_d_jac_smoother.o: amg_d_base_smoother_mod.o
amg_d_jac_smoother.o: amg_d_diag_solver.o
amg_dprecinit.o amg_dprecset.o: amg_d_diag_solver.o amg_d_ilu_solver.o \
		amg_d_umf_solver.o amg_d_as_smoother.o amg_d_jac_smoother.o \
	        amg_d_id_solver.o amg_d_slu_solver.o amg_d_sludist_solver.o

amg_s_mumps_solver.o amg_s_gs_solver.o amg_s_id_solver.o amg_s_slu_solver.o \
amg_s_diag_solver.o amg_s_ilu_solver.o: amg_s_base_solver_mod.o  amg_s_prec_type.o
amg_s_ilu_fact_mod.o: amg_base_prec_type.o  amg_s_base_solver_mod.o
amg_s_ilu_solver.o amg_s_iluk_fact.o: amg_s_ilu_fact_mod.o
amg_s_as_smoother.o amg_s_jac_smoother.o: amg_s_base_smoother_mod.o
amg_s_jac_smoother.o: amg_s_diag_solver.o
amg_sprecinit.o amg_sprecset.o: amg_s_diag_solver.o amg_s_ilu_solver.o \
		amg_s_as_smoother.o amg_s_jac_smoother.o \
	        amg_s_id_solver.o amg_s_slu_solver.o

amg_z_mumps_solver.o amg_z_gs_solver.o amg_z_id_solver.o amg_z_sludist_solver.o amg_z_slu_solver.o \
amg_z_umf_solver.o amg_z_diag_solver.o amg_z_ilu_solver.o: amg_z_base_solver_mod.o  amg_z_prec_type.o
amg_z_ilu_fact_mod.o: amg_base_prec_type.o  amg_z_base_solver_mod.o
amg_z_ilu_solver.o amg_z_iluk_fact.o: amg_z_ilu_fact_mod.o
amg_z_as_smoother.o amg_z_jac_smoother.o: amg_z_base_smoother_mod.o
amg_z_jac_smoother.o: amg_z_diag_solver.o
amg_zprecinit.o amg_zprecset.o: amg_z_diag_solver.o amg_z_ilu_solver.o \
		amg_z_umf_solver.o amg_z_as_smoother.o amg_z_jac_smoother.o \
	        amg_z_id_solver.o amg_z_slu_solver.o amg_z_sludist_solver.o

amg_c_mumps_solver.o amg_c_gs_solver.o amg_c_id_solver.o amg_c_sludist_solver.o amg_c_slu_solver.o \
amg_c_diag_solver.o amg_c_ilu_solver.o: amg_c_base_solver_mod.o  amg_c_prec_type.o
amg_c_ilu_fact_mod.o: amg_base_prec_type.o  amg_c_base_solver_mod.o
amg_c_ilu_solver.o amg_c_iluk_fact.o: amg_c_ilu_fact_mod.o
amg_c_as_smoother.o amg_c_jac_smoother.o: amg_c_base_smoother_mod.o
amg_c_jac_smoother.o: amg_c_diag_solver.o
amg_cprecinit.o amg_cprecset.o: amg_c_diag_solver.o amg_c_ilu_solver.o \
		amg_c_as_smoother.o amg_c_jac_smoother.o \
	        amg_c_id_solver.o amg_c_slu_solver.o amg_c_sludist_solver.o

amg_base_ainv_mod.o: amg_base_prec_type.o
amg_s_base_ainv_mod.o amg_d_base_ainv_mod.o amg_c_base_ainv_mod.o amg_z_base_ainv_mod.o: amg_base_ainv_mod.o
amg_s_ainv_solver.o: amg_base_ainv_mod.o amg_s_base_ainv_mod.o
amg_d_ainv_solver.o: amg_base_ainv_mod.o amg_d_base_ainv_mod.o
amg_c_ainv_solver.o: amg_base_ainv_mod.o amg_c_base_ainv_mod.o
amg_z_ainv_solver.o: amg_base_ainv_mod.o amg_z_base_ainv_mod.o
amg_s_invk_solver.o: amg_base_ainv_mod.o amg_s_base_ainv_mod.o
amg_d_invk_solver.o: amg_base_ainv_mod.o amg_d_base_ainv_mod.o
amg_c_invk_solver.o: amg_base_ainv_mod.o amg_c_base_ainv_mod.o
amg_z_invk_solver.o: amg_base_ainv_mod.o amg_z_base_ainv_mod.o
amg_s_invt_solver.o: amg_base_ainv_mod.o amg_s_base_ainv_mod.o
amg_d_invt_solver.o: amg_base_ainv_mod.o amg_d_base_ainv_mod.o
amg_c_invt_solver.o: amg_base_ainv_mod.o amg_c_base_ainv_mod.o
amg_z_invt_solver.o: amg_base_ainv_mod.o amg_z_base_ainv_mod.o
amg_ainv_mod.o: amg_s_ainv_solver.o amg_d_ainv_solver.o amg_c_ainv_solver.o \
	amg_z_ainv_solver.o amg_s_invk_solver.o amg_d_invk_solver.o amg_c_invk_solver.o \
	amg_z_invk_solver.o amg_s_invt_solver.o amg_d_invt_solver.o amg_c_invt_solver.o \
	amg_z_invt_solver.o

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: implclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS) *$(.mod)

implclean:
	cd impl && $(MAKE) clean

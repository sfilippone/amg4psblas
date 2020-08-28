include ../Make.inc
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
HERE=.

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)


DMODOBJS=mld_d_prec_type.o \
	mld_d_inner_mod.o mld_d_ilu_solver.o mld_d_diag_solver.o mld_d_jac_smoother.o mld_d_as_smoother.o \
	mld_d_umf_solver.o mld_d_slu_solver.o mld_d_sludist_solver.o mld_d_id_solver.o\
	mld_d_base_solver_mod.o mld_d_base_smoother_mod.o mld_d_onelev_mod.o \
	mld_d_gs_solver.o mld_d_mumps_solver.o \
	mld_d_base_aggregator_mod.o \
	mld_d_dec_aggregator_mod.o mld_d_symdec_aggregator_mod.o
#mld_d_bcmatch_aggregator_mod.o

SMODOBJS=mld_s_prec_type.o mld_s_ilu_fact_mod.o \
	mld_s_inner_mod.o mld_s_ilu_solver.o mld_s_diag_solver.o mld_s_jac_smoother.o mld_s_as_smoother.o \
	mld_s_slu_solver.o mld_s_id_solver.o\
	mld_s_base_solver_mod.o mld_s_base_smoother_mod.o mld_s_onelev_mod.o  \
	mld_s_gs_solver.o mld_s_mumps_solver.o \
	mld_s_base_aggregator_mod.o  \
	mld_s_dec_aggregator_mod.o mld_s_symdec_aggregator_mod.o

ZMODOBJS=mld_z_prec_type.o mld_z_ilu_fact_mod.o \
	mld_z_inner_mod.o mld_z_ilu_solver.o mld_z_diag_solver.o mld_z_jac_smoother.o mld_z_as_smoother.o \
	mld_z_umf_solver.o mld_z_slu_solver.o mld_z_sludist_solver.o mld_z_id_solver.o\
	mld_z_base_solver_mod.o mld_z_base_smoother_mod.o mld_z_onelev_mod.o \
	mld_z_gs_solver.o mld_z_mumps_solver.o \
	mld_z_base_aggregator_mod.o  \
	mld_z_dec_aggregator_mod.o mld_z_symdec_aggregator_mod.o

CMODOBJS=mld_c_prec_type.o mld_c_ilu_fact_mod.o \
	mld_c_inner_mod.o mld_c_ilu_solver.o mld_c_diag_solver.o mld_c_jac_smoother.o mld_c_as_smoother.o \
	mld_c_slu_solver.o mld_c_id_solver.o\
	mld_c_base_solver_mod.o mld_c_base_smoother_mod.o mld_c_onelev_mod.o  \
	mld_c_gs_solver.o mld_c_mumps_solver.o \
	mld_c_base_aggregator_mod.o  \
	mld_c_dec_aggregator_mod.o mld_c_symdec_aggregator_mod.o



MODOBJS=mld_base_prec_type.o mld_prec_type.o  mld_prec_mod.o \
	 mld_s_prec_mod.o  mld_d_prec_mod.o mld_c_prec_mod.o mld_z_prec_mod.o \
	$(SMODOBJS) $(DMODOBJS) $(CMODOBJS) $(ZMODOBJS) 


OBJS=$(MODOBJS) 

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=libmld_prec.a

all: lib impld

impld: $(OBJS)
	$(MAKE) -C impl

lib: $(OBJS) impld
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p mld_const.h $(INCDIR)
	/bin/cp -p *$(.mod) $(MODDIR)


$(MODOBJS): $(PSBLAS_MODDIR)/$(BASEMODNAME)$(.mod)

mld_base_prec_type.o: mld_const.h 
mld_s_prec_type.o mld_d_prec_type.o mld_c_prec_type.o mld_z_prec_type.o	: mld_base_prec_type.o
mld_prec_type.o: mld_s_prec_type.o mld_d_prec_type.o mld_c_prec_type.o mld_z_prec_type.o	
mld_prec_mod.o: mld_prec_type.o mld_s_prec_mod.o mld_d_prec_mod.o mld_c_prec_mod.o mld_z_prec_mod.o

$(SINNEROBJS) $(SOUTEROBJS): $(SMODOBJS)
$(DINNEROBJS) $(DOUTEROBJS): $(DMODOBJS)
$(CINNEROBJS) $(COUTEROBJS): $(CMODOBJS)
$(ZINNEROBJS) $(ZOUTEROBJS): $(ZMODOBJS)

mld_s_inner_mod.o: mld_s_prec_type.o
mld_d_inner_mod.o: mld_d_prec_type.o
mld_c_inner_mod.o: mld_c_prec_type.o
mld_z_inner_mod.o: mld_z_prec_type.o

mld_s_prec_mod.o: $(SMODOBJS) 
mld_d_prec_mod.o: $(DMODOBJS) 
mld_c_prec_mod.o: $(CMODOBJS) 
mld_z_prec_mod.o: $(ZMODOBJS) 


mld_s_prec_type.o: mld_s_onelev_mod.o
mld_d_prec_type.o: mld_d_onelev_mod.o
mld_c_prec_type.o: mld_c_onelev_mod.o
mld_z_prec_type.o: mld_z_onelev_mod.o

mld_s_onelev_mod.o: mld_s_base_smoother_mod.o mld_s_dec_aggregator_mod.o
mld_d_onelev_mod.o: mld_d_base_smoother_mod.o mld_d_dec_aggregator_mod.o
mld_c_onelev_mod.o: mld_c_base_smoother_mod.o mld_c_dec_aggregator_mod.o
mld_z_onelev_mod.o: mld_z_base_smoother_mod.o mld_z_dec_aggregator_mod.o

mld_s_base_aggregator_mod.o: mld_base_prec_type.o
mld_s_dec_aggregator_mod.o: mld_s_base_aggregator_mod.o
mld_s_hybrid_aggregator_mod.o mld_s_symdec_aggregator_mod.o: mld_s_dec_aggregator_mod.o

mld_d_base_aggregator_mod.o: mld_base_prec_type.o
mld_d_dec_aggregator_mod.o: mld_d_base_aggregator_mod.o
mld_d_hybrid_aggregator_mod.o mld_d_symdec_aggregator_mod.o: mld_d_dec_aggregator_mod.o

mld_c_base_aggregator_mod.o: mld_base_prec_type.o
mld_c_dec_aggregator_mod.o: mld_c_base_aggregator_mod.o
mld_c_hybrid_aggregator_mod.o mld_c_symdec_aggregator_mod.o: mld_c_dec_aggregator_mod.o

mld_z_base_aggregator_mod.o: mld_base_prec_type.o
mld_z_dec_aggregator_mod.o: mld_z_base_aggregator_mod.o
mld_z_hybrid_aggregator_mod.o mld_z_symdec_aggregator_mod.o: mld_z_dec_aggregator_mod.o

mld_s_base_smoother_mod.o: mld_s_base_solver_mod.o
mld_d_base_smoother_mod.o: mld_d_base_solver_mod.o
mld_c_base_smoother_mod.o: mld_c_base_solver_mod.o
mld_z_base_smoother_mod.o: mld_z_base_solver_mod.o


mld_s_base_solver_mod.o mld_d_base_solver_mod.o mld_c_base_solver_mod.o mld_z_base_solver_mod.o: mld_base_prec_type.o

mld_d_mumps_solver.o mld_d_gs_solver.o mld_d_id_solver.o mld_d_sludist_solver.o mld_d_slu_solver.o \
mld_d_umf_solver.o mld_d_diag_solver.o mld_d_ilu_solver.o: mld_d_base_solver_mod.o  mld_d_prec_type.o

#mld_d_ilu_fact_mod.o: mld_base_prec_type.o  mld_d_base_solver_mod.o
#mld_d_ilu_solver.o mld_d_iluk_fact.o: mld_d_ilu_fact_mod.o
mld_d_as_smoother.o mld_d_jac_smoother.o: mld_d_base_smoother_mod.o
mld_d_jac_smoother.o: mld_d_diag_solver.o
mld_dprecinit.o mld_dprecset.o: mld_d_diag_solver.o mld_d_ilu_solver.o \
		mld_d_umf_solver.o mld_d_as_smoother.o mld_d_jac_smoother.o \
	        mld_d_id_solver.o mld_d_slu_solver.o mld_d_sludist_solver.o

mld_s_mumps_solver.o mld_s_gs_solver.o mld_s_id_solver.o mld_s_slu_solver.o \
mld_s_diag_solver.o mld_s_ilu_solver.o: mld_s_base_solver_mod.o  mld_s_prec_type.o
mld_s_ilu_fact_mod.o: mld_base_prec_type.o  mld_s_base_solver_mod.o
mld_s_ilu_solver.o mld_s_iluk_fact.o: mld_s_ilu_fact_mod.o
mld_s_as_smoother.o mld_s_jac_smoother.o: mld_s_base_smoother_mod.o
mld_s_jac_smoother.o: mld_s_diag_solver.o
mld_sprecinit.o mld_sprecset.o: mld_s_diag_solver.o mld_s_ilu_solver.o \
		mld_s_as_smoother.o mld_s_jac_smoother.o \
	        mld_s_id_solver.o mld_s_slu_solver.o 

mld_z_mumps_solver.o mld_z_gs_solver.o mld_z_id_solver.o mld_z_sludist_solver.o mld_z_slu_solver.o \
mld_z_umf_solver.o mld_z_diag_solver.o mld_z_ilu_solver.o: mld_z_base_solver_mod.o  mld_z_prec_type.o
mld_z_ilu_fact_mod.o: mld_base_prec_type.o  mld_z_base_solver_mod.o	
mld_z_ilu_solver.o mld_z_iluk_fact.o: mld_z_ilu_fact_mod.o
mld_z_as_smoother.o mld_z_jac_smoother.o: mld_z_base_smoother_mod.o
mld_z_jac_smoother.o: mld_z_diag_solver.o
mld_zprecinit.o mld_zprecset.o: mld_z_diag_solver.o mld_z_ilu_solver.o \
		mld_z_umf_solver.o mld_z_as_smoother.o mld_z_jac_smoother.o \
	        mld_z_id_solver.o mld_z_slu_solver.o mld_z_sludist_solver.o

mld_c_mumps_solver.o mld_c_gs_solver.o mld_c_id_solver.o mld_c_sludist_solver.o mld_c_slu_solver.o \
mld_c_diag_solver.o mld_c_ilu_solver.o: mld_c_base_solver_mod.o  mld_c_prec_type.o
mld_c_ilu_fact_mod.o: mld_base_prec_type.o  mld_c_base_solver_mod.o
mld_c_ilu_solver.o mld_c_iluk_fact.o: mld_c_ilu_fact_mod.o
mld_c_as_smoother.o mld_c_jac_smoother.o: mld_c_base_smoother_mod.o
mld_c_jac_smoother.o: mld_c_diag_solver.o
mld_cprecinit.o mld_cprecset.o: mld_c_diag_solver.o mld_c_ilu_solver.o \
		mld_c_as_smoother.o mld_c_jac_smoother.o \
	        mld_c_id_solver.o mld_c_slu_solver.o mld_c_sludist_solver.o



veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: implclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS) *$(.mod)

implclean:
	$(MAKE) -C impl clean

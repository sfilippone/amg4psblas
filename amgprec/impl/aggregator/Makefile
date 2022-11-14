include ../../../Make.inc
LIBDIR=../../../lib
INCDIR=../../../include
MODDIR=../../../modules
HERE=../..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)
CXXINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(FMFLAG)/.

#CINCLUDES= -I${SUPERLU_INCDIR} -I${HSL_INCDIR}  -I${SPRAL_INCDIR} -I/home/users/pasqua/Ambra/BootCMatch/include -lBCM -L/home/users/pasqua/Ambra/BootCMatch/lib  -lm

FOBJS= \
amg_s_dec_aggregator_mat_asb.o \
amg_s_dec_aggregator_mat_bld.o \
amg_s_dec_aggregator_tprol.o \
amg_s_symdec_aggregator_tprol.o  \
amg_s_map_to_tprol.o amg_s_soc1_map_bld.o amg_s_soc2_map_bld.o\
amg_s_rap.o amg_s_ptap_bld.o \
amg_saggrmat_minnrg_bld.o\
amg_saggrmat_nosmth_bld.o amg_saggrmat_smth_bld.o \
amg_d_dec_aggregator_mat_asb.o \
amg_d_dec_aggregator_mat_bld.o \
amg_d_dec_aggregator_tprol.o \
amg_d_symdec_aggregator_tprol.o  \
amg_d_map_to_tprol.o amg_d_soc1_map_bld.o amg_d_soc2_map_bld.o \
amg_d_rap.o amg_d_ptap_bld.o \
amg_daggrmat_minnrg_bld.o \
amg_daggrmat_nosmth_bld.o amg_daggrmat_smth_bld.o \
amg_c_dec_aggregator_mat_asb.o \
amg_c_dec_aggregator_mat_bld.o \
amg_c_dec_aggregator_tprol.o \
amg_c_symdec_aggregator_tprol.o  \
amg_c_map_to_tprol.o amg_c_soc1_map_bld.o amg_c_soc2_map_bld.o\
amg_c_rap.o amg_c_ptap_bld.o \
amg_caggrmat_minnrg_bld.o\
amg_caggrmat_nosmth_bld.o amg_caggrmat_smth_bld.o \
amg_z_dec_aggregator_mat_asb.o \
amg_z_dec_aggregator_mat_bld.o \
amg_z_dec_aggregator_tprol.o \
amg_z_symdec_aggregator_tprol.o  \
amg_z_map_to_tprol.o amg_z_soc1_map_bld.o amg_z_soc2_map_bld.o\
amg_z_rap.o amg_z_ptap_bld.o \
amg_zaggrmat_minnrg_bld.o\
amg_zaggrmat_nosmth_bld.o amg_zaggrmat_smth_bld.o \
amg_d_parmatch_aggregator_mat_bld.o \
amg_d_parmatch_aggregator_mat_asb.o \
amg_d_parmatch_aggregator_inner_mat_asb.o \
amg_d_parmatch_aggregator_tprol.o \
amg_d_parmatch_spmm_bld.o \
amg_d_parmatch_spmm_bld_ov.o \
amg_d_parmatch_unsmth_bld.o \
amg_d_parmatch_smth_bld.o \
amg_d_parmatch_spmm_bld_inner.o \
amg_s_parmatch_aggregator_mat_bld.o \
amg_s_parmatch_aggregator_mat_asb.o \
amg_s_parmatch_aggregator_inner_mat_asb.o \
amg_s_parmatch_aggregator_tprol.o \
amg_s_parmatch_spmm_bld.o \
amg_s_parmatch_spmm_bld_ov.o \
amg_s_parmatch_unsmth_bld.o \
amg_s_parmatch_smth_bld.o \
amg_s_parmatch_spmm_bld_inner.o

MPCOBJS=MatchBoxPC.o \
sendBundledMessages.o \
initialize.o \
extractUChunk.o \
isAlreadyMatched.o \
findOwnerOfGhost.o \
computeCandidateMate.o \
parallelComputeCandidateMateB.o \
processMatchedVertices.o \
processMatchedVerticesAndSendMessages.o \
processCrossEdge.o \
queueTransfer.o \
processMessages.o \
processExposedVertex.o \
algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateC.o \
algoDistEdgeApproxDomEdgesLinearSearchMesgBndlSmallMateCMP.o

OBJS = $(FOBJS) $(MPCOBJS)

LIBNAME=libamg_prec.a

objs: $(OBJS)

lib:  objs
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

mpobjs:
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

!
!
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.7)
!
!    (C) Copyright 2020
!
!        Salvatore Filippone
!        Pasqua D'Ambra
!        Fabio Durastante
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
! File: amg_base_prec_type.F90
!
! Module: amg_base_prec_type
!
!  Constants and utilities in common to all type variants of MLD preconditioners.
!  - integer constants defining the preconditioner;
!  - character constants describing the preconditioner (used by the routines
!    printing out a preconditioner description);
!  - the interfaces to the routines for the management of the preconditioner
!    data structure (see below);
!  - The data type encapsulating the parameters defining the ML preconditioner
!  - The data type encapsulating the basic aggregation map.
!
!  It contains routines for
!  - converting character constants defining the preconditioner into integer
!    constants;
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.
!

module amg_base_prec_type

  !
  ! This reduces the size of .mod file. Without the ONLY clause compilation
  ! blows up on some systems.
  !
  use psb_const_mod
  use psb_base_mod, only :&
       & psb_desc_type, psb_ctxt_type,&
       & psb_ipk_, psb_dpk_, psb_spk_, psb_epk_,  &
       & psb_cdfree, psb_halo_, psb_none_, psb_sum_, psb_avg_, &
       & psb_nohalo_, psb_square_root_, psb_toupper, psb_root_,&
       & psb_sizeof_ip, psb_sizeof_lp, psb_sizeof_sp, &
       & psb_sizeof_dp, psb_sizeof,&
       & psb_cd_get_context, psb_info, psb_min, psb_sum, psb_bcast,&
       & psb_sizeof, psb_free, psb_cdfree, &
       & psb_errpush, psb_act_abort_, psb_act_ret_,&
       & psb_erractionsave, psb_erractionrestore, psb_error, psb_get_errstatus, &
       & psb_get_erraction, psb_success_, psb_err_alloc_dealloc_,&
       & psb_err_from_subroutine_, psb_err_missing_override_method_, &
       & psb_error_handler, psb_out_unit, psb_err_unit

  !
  ! Version numbers
  !
  character(len=*), parameter   :: amg_version_string_ = "1.1.0"
  integer(psb_ipk_), parameter  :: amg_version_major_  = 1
  integer(psb_ipk_), parameter  :: amg_version_minor_  = 1
  integer(psb_ipk_), parameter  :: amg_patchlevel_     = 0

  type amg_ml_parms
    integer(psb_ipk_) :: sweeps_pre, sweeps_post
    integer(psb_ipk_) :: ml_cycle
    integer(psb_ipk_) :: aggr_type, par_aggr_alg
    integer(psb_ipk_) :: aggr_ord, aggr_prol
    integer(psb_ipk_) :: aggr_omega_alg, aggr_eig, aggr_filter
    integer(psb_ipk_) :: coarse_mat, coarse_solve
  contains
    procedure, pass(pm) :: get_coarse_mat  => ml_parms_get_coarse_mat
    procedure, pass(pm) :: get_coarse  => ml_parms_get_coarse
    procedure, pass(pm) :: clone       => ml_parms_clone
    procedure, pass(pm) :: descr       => ml_parms_descr
    procedure, pass(pm) :: mlcycledsc  => ml_parms_mlcycledsc
    procedure, pass(pm) :: mldescr     => ml_parms_mldescr
    procedure, pass(pm) :: coarsedescr => ml_parms_coarsedescr
    procedure, pass(pm) :: printout    => ml_parms_printout
  end type amg_ml_parms


  type, extends(amg_ml_parms) :: amg_sml_parms
    real(psb_spk_) :: aggr_omega_val,  aggr_thresh
  contains
    procedure, pass(pm) :: clone => s_ml_parms_clone
    procedure, pass(pm) :: descr => s_ml_parms_descr
    procedure, pass(pm) :: printout => s_ml_parms_printout
  end type amg_sml_parms

  type, extends(amg_ml_parms) :: amg_dml_parms
    real(psb_dpk_) :: aggr_omega_val,  aggr_thresh
  contains
    procedure, pass(pm) :: clone => d_ml_parms_clone
    procedure, pass(pm) :: descr => d_ml_parms_descr
    procedure, pass(pm) :: printout => d_ml_parms_printout
  end type amg_dml_parms



  type amg_iaggr_data
    !
    ! Aggregation data and defaults:
    !
    ! 1. min_coarse_size = 0      Default target size will be computed  as
    !                               40*(N_fine)**(1./3.)
    !                             We are assuming that the coarse size fits in
    !                             integer range of psb_ipk_, but this is
    !                             not very restrictive
    integer(psb_ipk_)                  :: min_coarse_size = -ione
    integer(psb_ipk_)                  :: min_coarse_size_per_process = -ione
    integer(psb_lpk_)                  :: target_coarse_size
    ! 2. maximum number of levels.   Defaults to  20
    integer(psb_ipk_)                  :: max_levs    = 20_psb_ipk_
  contains
    procedure, pass(ag) :: default => i_ag_default
  end type amg_iaggr_data

  type, extends(amg_iaggr_data) :: amg_saggr_data
    ! 3. min_cr_ratio   = 1.5
    real(psb_spk_)                     :: min_cr_ratio   = 1.5_psb_spk_
    real(psb_spk_)                     :: op_complexity  = szero
    real(psb_spk_)                     :: avg_cr         = szero
  contains
    procedure, pass(ag) ::   default => s_ag_default
  end type amg_saggr_data

  type, extends(amg_iaggr_data) :: amg_daggr_data
    ! 3. min_cr_ratio   = 1.5
    real(psb_dpk_)                     :: min_cr_ratio   = 1.5_psb_dpk_
    real(psb_dpk_)                     :: op_complexity  = dzero
    real(psb_dpk_)                     :: avg_cr         = dzero
  contains
    procedure, pass(ag) ::   default => d_ag_default
  end type amg_daggr_data



  !
  ! Entries in iprcparm
  !
  ! These are in baseprec
  !
  integer(psb_ipk_), parameter :: amg_smoother_type_   =  1
  integer(psb_ipk_), parameter :: amg_sub_solve_       =  2
  integer(psb_ipk_), parameter :: amg_sub_restr_       =  3
  integer(psb_ipk_), parameter :: amg_sub_prol_        =  4
  integer(psb_ipk_), parameter :: amg_sub_ovr_         =  6
  integer(psb_ipk_), parameter :: amg_sub_fillin_      =  7
  integer(psb_ipk_), parameter :: amg_ilu_scale_       =  8

  !
  ! These are in onelev
  !
  integer(psb_ipk_), parameter :: amg_ml_cycle_             = 20
  integer(psb_ipk_), parameter :: amg_smoother_sweeps_pre_  = 21
  integer(psb_ipk_), parameter :: amg_smoother_sweeps_post_ = 22
  integer(psb_ipk_), parameter :: amg_aggr_type_            = 23
  integer(psb_ipk_), parameter :: amg_aggr_prol_            = 24
  integer(psb_ipk_), parameter :: amg_par_aggr_alg_         = 25
  integer(psb_ipk_), parameter :: amg_aggr_ord_             = 26
  integer(psb_ipk_), parameter :: amg_aggr_omega_alg_       = 27
  integer(psb_ipk_), parameter :: amg_aggr_eig_             = 28
  integer(psb_ipk_), parameter :: amg_aggr_filter_          = 29
  integer(psb_ipk_), parameter :: amg_coarse_mat_           = 30
  integer(psb_ipk_), parameter :: amg_coarse_solve_         = 31
  integer(psb_ipk_), parameter :: amg_coarse_sweeps_        = 32
  integer(psb_ipk_), parameter :: amg_coarse_fillin_        = 33
  integer(psb_ipk_), parameter :: amg_coarse_subsolve_      = 34
  integer(psb_ipk_), parameter :: amg_smoother_sweeps_      = 36
  integer(psb_ipk_), parameter :: amg_solver_sweeps_        = 37
  integer(psb_ipk_), parameter :: amg_min_coarse_size_      = 38
  integer(psb_ipk_), parameter :: amg_n_prec_levs_          = 39
  integer(psb_ipk_), parameter :: amg_max_levs_             = 40
  integer(psb_ipk_), parameter :: amg_min_cr_ratio_         = 41
  integer(psb_ipk_), parameter :: amg_outer_sweeps_         = 42
  integer(psb_ipk_), parameter :: amg_ifpsz_                = 43

  !
  ! Legal values for entry: amg_smoother_type_
  !
  integer(psb_ipk_), parameter :: amg_min_prec_ = 0
  integer(psb_ipk_), parameter :: amg_noprec_   = 0
  integer(psb_ipk_), parameter :: amg_base_smooth_ = 0
  integer(psb_ipk_), parameter :: amg_jac_      = 1
  integer(psb_ipk_), parameter :: amg_l1_jac_   = 2
  integer(psb_ipk_), parameter :: amg_bjac_     = 3
  integer(psb_ipk_), parameter :: amg_l1_bjac_  = 4
  integer(psb_ipk_), parameter :: amg_as_       = 5
  integer(psb_ipk_), parameter :: amg_fbgs_     = 6
  integer(psb_ipk_), parameter :: amg_l1_gs_    = 7
  integer(psb_ipk_), parameter :: amg_l1_fbgs_  = 8
  integer(psb_ipk_), parameter :: amg_poly_     = 9
  integer(psb_ipk_), parameter :: amg_max_prec_ = 9
  !
  ! Constants for pre/post signaling. Now only used internally
  !
  integer(psb_ipk_), parameter :: amg_smooth_pre_     = 1
  integer(psb_ipk_), parameter :: amg_smooth_post_    = 2
  integer(psb_ipk_), parameter :: amg_smooth_both_    = 3

  !
  !  This is a quick&dirty fix, but I have nothing better now...
  !
  ! Legal values for entry: amg_sub_solve_
  !
  integer(psb_ipk_), parameter :: amg_slv_delta_     = amg_max_prec_+1
  integer(psb_ipk_), parameter :: amg_f_none_        = amg_slv_delta_+0
  integer(psb_ipk_), parameter :: amg_diag_scale_    = amg_slv_delta_+1
  integer(psb_ipk_), parameter :: amg_l1_diag_scale_ = amg_slv_delta_+2
  integer(psb_ipk_), parameter :: amg_gs_            = amg_slv_delta_+3
  integer(psb_ipk_), parameter :: amg_ilu_n_         = amg_slv_delta_+4
  integer(psb_ipk_), parameter :: amg_milu_n_        = amg_slv_delta_+5
  integer(psb_ipk_), parameter :: amg_ilu_t_         = amg_slv_delta_+6
  integer(psb_ipk_), parameter :: amg_slu_           = amg_slv_delta_+7
  integer(psb_ipk_), parameter :: amg_umf_           = amg_slv_delta_+8
  integer(psb_ipk_), parameter :: amg_sludist_       = amg_slv_delta_+9
  integer(psb_ipk_), parameter :: amg_mumps_         = amg_slv_delta_+10
  integer(psb_ipk_), parameter :: amg_bwgs_          = amg_slv_delta_+11
  integer(psb_ipk_), parameter :: amg_krm_           = amg_slv_delta_+12
  integer(psb_ipk_), parameter :: amg_max_sub_solve_ = amg_slv_delta_+12
  integer(psb_ipk_), parameter :: amg_min_sub_solve_ = amg_diag_scale_

  !
  ! Legal values for entry: amg_ilu_scale_
  !
  integer(psb_ipk_), parameter :: amg_ilu_scale_none_    = 0
  integer(psb_ipk_), parameter :: amg_ilu_scale_maxval_  = 1
  integer(psb_ipk_), parameter :: amg_ilu_scale_diag_    = 2
  integer(psb_ipk_), parameter :: amg_ilu_scale_arwsum_  = 3
  integer(psb_ipk_), parameter :: amg_ilu_scale_aclsum_  = 4
  integer(psb_ipk_), parameter :: amg_ilu_scale_arcsum_  = 5
  ! For the time being enable only maxval scale
  integer(psb_ipk_), parameter :: amg_max_ilu_scale_     = 1
  !
  ! Legal values for entry: amg_ml_cycle_
  !
  integer(psb_ipk_), parameter :: amg_no_ml_        = 0
  integer(psb_ipk_), parameter :: amg_add_ml_       = 1
  integer(psb_ipk_), parameter :: amg_mult_ml_      = 2
  integer(psb_ipk_), parameter :: amg_vcycle_ml_    = 3
  integer(psb_ipk_), parameter :: amg_wcycle_ml_    = 4
  integer(psb_ipk_), parameter :: amg_kcycle_ml_    = 5
  integer(psb_ipk_), parameter :: amg_kcyclesym_ml_ = 6
  integer(psb_ipk_), parameter :: amg_new_ml_prec_  = 7
  integer(psb_ipk_), parameter :: amg_mult_dev_ml_  = 7
  integer(psb_ipk_), parameter :: amg_max_ml_cycle_  = 8
  !
  ! Legal values for entry: amg_par_aggr_alg_
  !
  integer(psb_ipk_), parameter :: amg_dec_aggr_      = 0
  integer(psb_ipk_), parameter :: amg_sym_dec_aggr_  = 1
  integer(psb_ipk_), parameter :: amg_ext_aggr_      = 2
  integer(psb_ipk_), parameter :: amg_coupled_aggr_  = 3
  integer(psb_ipk_), parameter :: amg_max_par_aggr_alg_  = amg_coupled_aggr_
  !
  ! Legal values for entry: amg_aggr_type_
  !
  integer(psb_ipk_), parameter :: amg_noalg_       = 0
  integer(psb_ipk_), parameter :: amg_soc1_        = 1
  integer(psb_ipk_), parameter :: amg_soc2_        = 2
  integer(psb_ipk_), parameter :: amg_matchboxp_   = 3
  !
  ! Legal values for entry: amg_aggr_prol_
  !
  integer(psb_ipk_), parameter :: amg_no_smooth_   = 0
  integer(psb_ipk_), parameter :: amg_smooth_prol_ = 1
  integer(psb_ipk_), parameter :: amg_min_energy_  = 2
  ! Disabling min_energy  for the time being.
  integer(psb_ipk_), parameter :: amg_max_aggr_prol_=amg_smooth_prol_
  !
  ! Legal values for entry: amg_aggr_filter_
  !
  integer(psb_ipk_), parameter :: amg_no_filter_mat_  = 0
  integer(psb_ipk_), parameter :: amg_filter_mat_     = 1
  integer(psb_ipk_), parameter :: amg_max_filter_mat_ = amg_filter_mat_
  !
  ! Legal values for entry: amg_aggr_ord_
  !
  integer(psb_ipk_), parameter :: amg_aggr_ord_nat_      = 0
  integer(psb_ipk_), parameter :: amg_aggr_ord_desc_deg_ = 1
  integer(psb_ipk_), parameter :: amg_max_aggr_ord_      = amg_aggr_ord_desc_deg_
  !
  ! Legal values for entry: amg_aggr_omega_alg_
  !
  integer(psb_ipk_), parameter :: amg_eig_est_     = 0
  integer(psb_ipk_), parameter :: amg_user_choice_ = 999
  !
  ! Legal values for entry: amg_aggr_eig_
  !
  integer(psb_ipk_), parameter :: amg_max_norm_ = 0
  !
  ! Legal values for entry: amg_coarse_mat_
  !
  integer(psb_ipk_), parameter :: amg_distr_mat_      = 0
  integer(psb_ipk_), parameter :: amg_repl_mat_       = 1
  integer(psb_ipk_), parameter :: amg_max_coarse_mat_ = amg_repl_mat_
  !
  ! Legal values for entry: amg_poly_variant_
  !
  integer(psb_ipk_), parameter :: amg_poly_lottes_      = 0
  integer(psb_ipk_), parameter :: amg_poly_lottes_beta_ = 1
  integer(psb_ipk_), parameter :: amg_poly_new_         = 2
  
  integer(psb_ipk_), parameter :: amg_poly_rho_est_power_ = 0

  !
  ! Legal values for entry: amg_prec_status_
  !
  integer(psb_ipk_), parameter :: amg_prec_built_ = 98765

  !
  ! Entries in rprcparm: ILU(k,t) threshold, smoothed aggregation omega
  !
  integer(psb_ipk_), parameter :: amg_sub_iluthrs_    = 1
  integer(psb_ipk_), parameter :: amg_aggr_omega_val_ = 2
  integer(psb_ipk_), parameter :: amg_aggr_thresh_    = 3
  integer(psb_ipk_), parameter :: amg_coarse_iluthrs_ = 4
  integer(psb_ipk_), parameter :: amg_solver_eps_     = 6
  integer(psb_ipk_), parameter :: amg_rfpsz_          = 8
  !
  ! Is the current solver local or global
  !
  integer(psb_ipk_), parameter :: amg_local_solver_   = 0
  integer(psb_ipk_), parameter :: amg_global_solver_  = 1

  !
  ! Entries for mumps
  !
  ! Size of the control vectors
  integer, parameter :: amg_mumps_icntl_size=40
  integer, parameter :: amg_mumps_rcntl_size=15

  !
  ! Fields for sparse matrices ensembles stored in av()
  !
  integer(psb_ipk_), parameter :: amg_l_pr_        = 1
  integer(psb_ipk_), parameter :: amg_u_pr_        = 2
  integer(psb_ipk_), parameter :: amg_bp_ilu_avsz_ = 2
  integer(psb_ipk_), parameter :: amg_ap_nd_       = 3
  integer(psb_ipk_), parameter :: amg_ac_          = 4
  integer(psb_ipk_), parameter :: amg_sm_pr_t_     = 5
  integer(psb_ipk_), parameter :: amg_sm_pr_       = 6
  integer(psb_ipk_), parameter :: amg_smth_avsz_   = 6
  integer(psb_ipk_), parameter :: amg_max_avsz_    = amg_smth_avsz_

  !
  ! Character constants used by amg_file_prec_descr
  !
  character(len=19), parameter, private :: &
       &  eigen_estimates(0:0)=(/'infinity norm     '/)
  character(len=15), parameter, private :: &
       &  aggr_prols(0:3)=(/'unsmoothed    ','smoothed      ',&
       &           'min energy    ','bizr. smoothed'/)
  character(len=15), parameter, private :: &
       &  aggr_filters(0:1)=(/'no filtering  ','filtering     '/)
  character(len=15), parameter, private :: &
       &  matrix_names(0:1)=(/'distributed   ','replicated    '/)
  character(len=18), parameter, private :: &
       &  aggr_type_names(0:3)=(/'None              ',&
       &  'SOC measure 1     ', 'SOC Measure 2     ',&
       &  'Parallel Matching '/)
  character(len=18), parameter, private :: &
       &  par_aggr_alg_names(0:3)=(/&
       & 'decoupled aggr.   ', 'sym. dec. aggr.   ',&
       & 'user defined aggr.', 'coupled aggr.     '/)
  character(len=18), parameter, private :: &
       &  ord_names(0:1)=(/'Natural ordering  ','Desc. degree ord. '/)
  character(len=6), parameter, private :: &
       &  restrict_names(0:4)=(/'none ','halo ','     ','     ','     '/)
  character(len=12), parameter, private :: &
       &  prolong_names(0:3)=(/'none       ','sum        ', &
       & 'average    ','square root'/)
  character(len=15), parameter, private :: &
       &  ml_names(0:7)=(/'none          ','additive      ',&
       &  'multiplicative', 'VCycle        ','WCycle        ',&
       &  'KCycle        ','KCycleSym     ','new ML        '/)
  character(len=16), parameter :: &
       &  amg_fact_names(0:amg_max_sub_solve_)=(/&
       & 'none          ','Jacobi        ',&
       & 'L1-Jacobi     ','none          ','none          ',&
       & 'none          ','none          ','L1-GS         ',&
       & 'L1-FBGS       ','Polynomial    ','none          ','Point Jacobi  ',&
       & 'L1-Jacobi     ','Gauss-Seidel  ','ILU(n)        ',&
       & 'MILU(n)       ','ILU(t,n)      ',&
       & 'SuperLU       ','UMFPACK LU    ',&
       & 'SuperLU_Dist  ','MUMPS         ',&
       & 'Backward GS   ','Krylov Method '/)

  interface amg_check_def
    module procedure amg_icheck_def, amg_scheck_def, amg_dcheck_def
  end interface

  interface psb_bcast
    module procedure amg_ml_bcast, amg_sml_bcast, amg_dml_bcast
  end interface psb_bcast

  interface amg_equal_aggregation
    module procedure amg_d_equal_aggregation, amg_s_equal_aggregation
  end interface amg_equal_aggregation

  !
  ! default to no remapping.
  ! Will need a more sophisticated strategy.
  !
  logical, private, save :: do_remap=.false.

contains

  function amg_get_do_remap() result(res)
    implicit none
    logical :: res

    res = do_remap
  end function amg_get_do_remap

  subroutine amg_set_do_remap(val)
    implicit none
    logical, intent(in) :: val

    do_remap = val
  end subroutine amg_set_do_remap

  !
  ! Function: amg_stringval
  !
  !  This routine converts the string contained into string into the corresponding
  !  integer value.
  !
  ! Arguments:
  !    string  -  character(len=*), input.
  !               The string to be converted.
  !    val     -  integer, output.
  !               The integer value corresponding to the string
  !
  function amg_stringval(string) result(val)
    use psb_prec_const_mod
    implicit none
  ! Arguments
    character(len=*), intent(in) :: string
    integer(psb_ipk_) :: val
    character(len=*), parameter :: name='amg_stringval'
  ! Local variable
    integer :: index_tab
    character(len=128) ::string2
    index_tab=index(string,char(9))
    if (index_tab.NE.0)  then
      string2=string(1:index_tab-1)
    else
      string2=string
    endif
    select case(psb_toupper(trim(string2)))
    case('NONE')
      val = 0
    case('HALO')
      val = psb_halo_
    case('SUM')
      val = psb_sum_
    case('AVG')
      val = psb_avg_
    case('FACT_NONE')
      val = amg_f_none_
    case('FBGS')
      val = amg_fbgs_
    case('GS','FGS','FWGS')
      val = amg_gs_
    case('BGS','BWGS')
      val = amg_bwgs_
    case('ILU')
      val = amg_ilu_n_
    case('MILU')
      val = amg_milu_n_
    case('ILUT')
      val = amg_ilu_t_
    case('MUMPS')
      val = amg_mumps_
    case('UMF')
      val = amg_umf_
    case('SLU')
      val = amg_slu_
    case('SLUDIST')
      val = amg_sludist_
    case('DIAG')
      val = amg_diag_scale_
    case('L1-DIAG')
      val = amg_l1_diag_scale_
    case('ADD')
      val = amg_add_ml_
    case('MULT_DEV')
      val = amg_mult_dev_ml_
    case('MULT')
      val = amg_mult_ml_
    case('VCYCLE')
      val = amg_vcycle_ml_
    case('WCYCLE')
      val = amg_wcycle_ml_
    case('KCYCLE')
      val = amg_kcycle_ml_
    case('KCYCLESYM')
      val = amg_kcyclesym_ml_
    case('SOC2')
      val = amg_soc2_
    case('SOC1')
      val = amg_soc1_
    case('MATCHBOXP','PARMATCH')
      val = amg_matchboxp_
    case('COUPLED','COUP')
      val = amg_coupled_aggr_
    case('DEC','DECOUPLED')
      val = amg_dec_aggr_
    case('SYMDEC')
      val = amg_sym_dec_aggr_
    case('NAT','NATURAL')
      val =  amg_aggr_ord_nat_
    case('DESC','RDEGREE','DEGREE')
      val = amg_aggr_ord_desc_deg_
    case('REPL')
      val = amg_repl_mat_
    case('DIST')
      val = amg_distr_mat_
    case('UNSMOOTHED','NONSMOOTHED')
      val = amg_no_smooth_
    case('SMOOTHED')
      val = amg_smooth_prol_
    case('MINENERGY')
      val = amg_min_energy_
    case('NOPREC')
      val = amg_noprec_
    case('BJAC')
      val = amg_bjac_
    case('L1-GS')
      val = amg_l1_gs_
    case('L1-FBGS')
      val = amg_l1_fbgs_
    case('L1-BJAC')
      val = amg_l1_bjac_
    case('JAC','JACOBI')
      val = amg_jac_
    case('L1-JACOBI')
      val = amg_l1_jac_
    case('KRM')
      val = amg_krm_
    case('AS')
      val = amg_as_
    case('POLY')
      val = amg_poly_
    case('POLY_LOTTES')
      val = amg_poly_lottes_
    case('POLY_LOTTES_BETA')
      val = amg_poly_lottes_beta_
    case('POLY_NEW')
      val = amg_poly_new_
    case('POLY_RHO_EST_POWER')
      val =  amg_poly_rho_est_power_
    case('A_NORMI')
      val = amg_max_norm_
    case('USER_CHOICE')
      val = amg_user_choice_
    case('EIG_EST')
      val = amg_eig_est_
    case('FILTER')
      val = amg_filter_mat_
    case('NOFILTER','NO_FILTER')
      val = amg_no_filter_mat_
    case('OUTER_SWEEPS')
      val = amg_outer_sweeps_
    case('LOCAL_SOLVER')
      val = amg_local_solver_
    case('GLOBAL_SOLVER')
      val = amg_global_solver_
    case default
      val  = -1
    end select
  end function amg_stringval

  function amg_get_coarse_mat_name(val) result(res)
    character(len=15) :: res
    integer :: val
    select case(val)
    case (0,1)
      res = matrix_names(val)
    case default
      res = 'Unknown '
    end select
  end function amg_get_coarse_mat_name

  subroutine amg_warn_coarse_mat(val,expected)
    integer(psb_ipk_) :: val, expected
    if (val /= expected) then
      write(0,*) 'Warning: resetting COARSE_MAT on an existing hierarchy from ',&
           & amg_get_coarse_mat_name(val), ' to ',amg_get_coarse_mat_name(expected)
    end if
  end subroutine amg_warn_coarse_mat

  
  function  ml_parms_get_coarse_mat(pm) result(res)
    implicit none
    class(amg_ml_parms), intent(in)    :: pm
    integer(psb_ipk_) :: res
    res   = pm%coarse_mat
  end function ml_parms_get_coarse_mat

  subroutine  ml_parms_get_coarse(pm,pmin)
    implicit none
    class(amg_ml_parms), intent(inout) :: pm
    class(amg_ml_parms), intent(in)    :: pmin
    pm%coarse_mat   = pmin%coarse_mat
    pm%coarse_solve = pmin%coarse_solve
  end subroutine ml_parms_get_coarse



  subroutine ml_parms_printout(pm,iout)
    implicit none
    class(amg_ml_parms), intent(in) :: pm
    integer(psb_ipk_), intent(in)             :: iout

    write(iout,*) 'ML    : ',pm%ml_cycle
    write(iout,*) 'Sweeps: ',pm%sweeps_pre,pm%sweeps_post
    write(iout,*) 'AGGR  : ',pm%par_aggr_alg,pm%aggr_prol, pm%aggr_ord
    write(iout,*) '      : ',pm%aggr_omega_alg,pm%aggr_eig,pm%aggr_filter
    write(iout,*) 'COARSE: ',pm%coarse_mat,pm%coarse_solve
  end subroutine ml_parms_printout


  subroutine s_ml_parms_printout(pm,iout)
    implicit none
    class(amg_sml_parms), intent(in) :: pm
    integer(psb_ipk_), intent(in)             :: iout

    call pm%amg_ml_parms%printout(iout)
    write(iout,*) 'REAL  : ',pm%aggr_omega_val,pm%aggr_thresh
  end subroutine s_ml_parms_printout


  subroutine d_ml_parms_printout(pm,iout)
    implicit none
    class(amg_dml_parms), intent(in) :: pm
    integer(psb_ipk_), intent(in)             :: iout

    call pm%amg_ml_parms%printout(iout)
    write(iout,*) 'REAL  : ',pm%aggr_omega_val,pm%aggr_thresh
  end subroutine d_ml_parms_printout


  !
  ! Routines printing out a description of the preconditioner
  !
  subroutine ml_parms_mlcycledsc(pm,iout,info)

    Implicit None

    ! Arguments
    class(amg_ml_parms), intent(in) :: pm
    integer(psb_ipk_), intent(in)             :: iout
    integer(psb_ipk_), intent(out)            :: info
    info = psb_success_
    if ((pm%ml_cycle>=amg_no_ml_).and.(pm%ml_cycle<=amg_max_ml_cycle_)) then


      write(iout,*) '  Multilevel cycle: ',&
           &   ml_names(pm%ml_cycle)
      select case (pm%ml_cycle)
      case (amg_add_ml_)
        write(iout,*) '  Number of smoother sweeps : ',&
             & pm%sweeps_pre
      case (amg_mult_ml_,amg_vcycle_ml_, amg_wcycle_ml_, amg_kcycle_ml_, amg_kcyclesym_ml_)
        write(iout,*) '  Number of smoother sweeps : pre: ',&
             &  pm%sweeps_pre ,'  post: ', pm%sweeps_post
      end select

    end if
  end subroutine ml_parms_mlcycledsc

  subroutine ml_parms_mldescr(pm,iout,info,prefix)

    Implicit None

    ! Arguments
    class(amg_ml_parms), intent(in)        :: pm
    integer(psb_ipk_), intent(in)          :: iout
    integer(psb_ipk_), intent(out)         :: info
    character(len=*), intent(in), optional :: prefix

    character(1024)  :: prefix_ 
    info = psb_success_
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if
    
    if ((pm%ml_cycle>=amg_no_ml_).and.(pm%ml_cycle<=amg_max_ml_cycle_)) then


      write(iout,*) trim(prefix),'  Parallel aggregation algorithm: ',&
           &   par_aggr_alg_names(pm%par_aggr_alg)
      if (pm%aggr_type>0) write(iout,*) trim(prefix),'  Aggregation type: ',&
           & aggr_type_names(pm%aggr_type)
      !if (pm%par_aggr_alg /= amg_ext_aggr_) then
        if ( pm%aggr_ord /= amg_aggr_ord_nat_) &
             & write(iout,*) trim(prefix),'               with initial ordering: ',&
             &   ord_names(pm%aggr_ord)
        write(iout,*) trim(prefix),'  Aggregation prolongator: ', &
             &  aggr_prols(pm%aggr_prol)
        if (pm%aggr_prol /= amg_no_smooth_) then
          write(iout,*) trim(prefix),'              with: ', aggr_filters(pm%aggr_filter)
          if (pm%aggr_omega_alg == amg_eig_est_) then
            write(iout,*) trim(prefix),'  Damping omega computation: spectral radius estimate'
            write(iout,*) trim(prefix),'  Spectral radius estimate: ', &
                 & eigen_estimates(pm%aggr_eig)
          else if (pm%aggr_omega_alg == amg_user_choice_) then
            write(iout,*) trim(prefix),'  Damping omega computation: user defined value.'
          else
            write(iout,*) trim(prefix),'  Damping omega computation: unknown value in iprcparm!!'
          end if
        end if
      !end if
    else
      write(iout,*) trim(prefix),'  Multilevel type: Unkonwn value. Something is amiss....',&
           & pm%ml_cycle
    end if

    return

  end subroutine ml_parms_mldescr

  subroutine ml_parms_descr(pm,iout,info,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_ml_parms), intent(in)        :: pm
    integer(psb_ipk_), intent(in)          :: iout
    integer(psb_ipk_), intent(out)         :: info
    logical, intent(in), optional          :: coarse
    character(len=*), intent(in), optional :: prefix
    logical :: coarse_

    info = psb_success_
    if (present(coarse)) then
      coarse_ = coarse
    else
      coarse_ = .false.
    end if

    if (coarse_) then
      call pm%coarsedescr(iout,info,prefix=prefix)
    end if

    return

  end subroutine ml_parms_descr


  subroutine ml_parms_coarsedescr(pm,iout,info,prefix)


    Implicit None

    ! Arguments
    class(amg_ml_parms), intent(in)       :: pm
    integer(psb_ipk_), intent(in)         :: iout
    integer(psb_ipk_), intent(out)        :: info
    character(len=*), intent(in), optional :: prefix

    character(1024)  :: prefix_ 

    info = psb_success_
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if

    write(iout,*) trim(prefix),'  Coarse matrix: ',&
         & matrix_names(pm%coarse_mat)
    select case(pm%coarse_solve)
    case (amg_bjac_,amg_as_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'Block Jacobi'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case (amg_l1_bjac_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'L1-Block Jacobi'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case (amg_jac_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'Point Jacobi'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case (amg_l1_jac_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'L1-Jacobi'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case (amg_l1_fbgs_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'L1 Forward-Backward Gauss-Seidel (Hybrid)'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case (amg_l1_gs_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'L1 Gauss-Seidel (Hybrid)'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case (amg_fbgs_)
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & 'Forward-Backward Gauss-Seidel (Hybrid)'
      write(iout,*) trim(prefix),'  Number of sweeps : ',&
           & pm%sweeps_pre
    case default
      write(iout,*) trim(prefix),'  Coarse solver: ',&
           & amg_fact_names(pm%coarse_solve)
    end select
    
  end subroutine ml_parms_coarsedescr

  subroutine s_ml_parms_descr(pm,iout,info,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_sml_parms), intent(in)       :: pm
    integer(psb_ipk_), intent(in)          :: iout
    integer(psb_ipk_), intent(out)         :: info
    logical, intent(in), optional          :: coarse
    character(len=*), intent(in), optional :: prefix

    character(1024)  :: prefix_ 

    info = psb_success_
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if

    call pm%amg_ml_parms%descr(iout,info,coarse,prefix=prefix)
    if (pm%aggr_prol /= amg_no_smooth_) then
      write(iout,*) trim(prefix),'  Damping omega value  :',pm%aggr_omega_val
    end if
    write(iout,*) trim(prefix),'  Aggregation threshold:',pm%aggr_thresh

    return

  end subroutine s_ml_parms_descr

  subroutine d_ml_parms_descr(pm,iout,info,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_dml_parms), intent(in)       :: pm
    integer(psb_ipk_), intent(in)          :: iout
    integer(psb_ipk_), intent(out)         :: info
    logical, intent(in), optional          :: coarse
    character(len=*), intent(in), optional :: prefix

    character(1024)  :: prefix_ 

    info = psb_success_
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if

    call pm%amg_ml_parms%descr(iout,info,coarse,prefix=prefix)
    if (pm%aggr_prol /= amg_no_smooth_) then
      write(iout,*) trim(prefix),'  Damping omega value  :',pm%aggr_omega_val
    end if
    write(iout,*) trim(prefix),'  Aggregation threshold:',pm%aggr_thresh

    return

  end subroutine d_ml_parms_descr


  !
  ! Functions/subroutines checking if the preconditioner is correctly defined
  !

  function is_legal_base_prec(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_base_prec

    is_legal_base_prec = ((ip>=amg_noprec_).and.(ip<=amg_max_prec_))
    return
  end function is_legal_base_prec
  function is_int_non_negative(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_int_non_negative

    is_int_non_negative = (ip >= 0)
    return
  end function is_int_non_negative
  function is_legal_ilu_scale(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ilu_scale
    is_legal_ilu_scale = ((ip >= amg_ilu_scale_none_).and.(ip <= amg_max_ilu_scale_))
    return
  end function is_legal_ilu_scale
  function is_int_positive(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_int_positive

    is_int_positive = (ip >= 1)
    return
  end function is_int_positive
  function is_legal_prolong(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_prolong
    is_legal_prolong = ((ip>=psb_none_).and.(ip<=psb_square_root_))
    return
  end function is_legal_prolong
  function is_legal_restrict(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_restrict
    is_legal_restrict = ((ip == psb_nohalo_).or.(ip==psb_halo_))
    return
  end function is_legal_restrict
  function is_legal_ml_cycle(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_cycle

    is_legal_ml_cycle = ((ip>=amg_no_ml_).and.(ip<=amg_max_ml_cycle_))
    return
  end function is_legal_ml_cycle
  function is_legal_coupled_par_aggr_alg(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_coupled_par_aggr_alg

    is_legal_coupled_par_aggr_alg = (ip == amg_coupled_aggr_)
    return
  end function is_legal_coupled_par_aggr_alg
  function is_legal_decoupled_par_aggr_alg(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_decoupled_par_aggr_alg

    is_legal_decoupled_par_aggr_alg = ((ip>=amg_dec_aggr_).and.(ip<=amg_max_par_aggr_alg_))
    return
  end function is_legal_decoupled_par_aggr_alg
  function is_legal_ml_aggr_type(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_aggr_type

    is_legal_ml_aggr_type = (ip >= amg_soc1_) .and.  (ip <= amg_soc2_)
    return
  end function is_legal_ml_aggr_type
  function is_legal_ml_aggr_ord(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_aggr_ord

    is_legal_ml_aggr_ord = ((amg_aggr_ord_nat_<=ip).and.(ip<=amg_max_aggr_ord_))
    return
  end function is_legal_ml_aggr_ord
  function is_legal_ml_aggr_omega_alg(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_aggr_omega_alg

    is_legal_ml_aggr_omega_alg = ((ip == amg_eig_est_).or.(ip==amg_user_choice_))
    return
  end function is_legal_ml_aggr_omega_alg
  function is_legal_ml_aggr_eig(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_aggr_eig

    is_legal_ml_aggr_eig = (ip == amg_max_norm_)
    return
  end function is_legal_ml_aggr_eig
  function is_legal_ml_aggr_prol(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_aggr_prol

    is_legal_ml_aggr_prol = ((ip>=0).and.(ip<=amg_max_aggr_prol_))
    return
  end function is_legal_ml_aggr_prol
  function is_legal_ml_coarse_mat(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_coarse_mat

    is_legal_ml_coarse_mat = ((ip>=0).and.(ip<=amg_max_coarse_mat_))
    return
  end function is_legal_ml_coarse_mat
  function is_legal_aggr_filter(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_aggr_filter

    is_legal_aggr_filter = ((ip>=0).and.(ip<=amg_max_filter_mat_))
    return
  end function is_legal_aggr_filter
  function is_distr_ml_coarse_mat(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_distr_ml_coarse_mat

    is_distr_ml_coarse_mat = (ip == amg_distr_mat_)
    return
  end function is_distr_ml_coarse_mat
  function is_legal_ml_fact(ip)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ml_fact
    ! Here the minimum is really 1, amg_fact_none_ is not acceptable.
    is_legal_ml_fact = ((ip>=amg_min_sub_solve_)&
         & .and.(ip<=amg_max_sub_solve_))
    return
  end function is_legal_ml_fact
  function is_legal_ilu_fact(ip)
    use psb_prec_const_mod
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: is_legal_ilu_fact

    is_legal_ilu_fact = ((ip==amg_ilu_n_).or.&
         & (ip==amg_milu_n_).or.(ip==amg_ilu_t_))
    return
  end function is_legal_ilu_fact
  function is_legal_d_omega(ip)
    implicit none
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_d_omega
    is_legal_d_omega = ((ip>=0.0d0).and.(ip<=2.0d0))
    return
  end function is_legal_d_omega
  function is_legal_d_fact_thrs(ip)
    implicit none
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_d_fact_thrs

    is_legal_d_fact_thrs = (ip>=0.0d0)
    return
  end function is_legal_d_fact_thrs
  function is_legal_d_aggr_thrs(ip)
    implicit none
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_d_aggr_thrs

    is_legal_d_aggr_thrs = (ip>=0.0d0)
    return
  end function is_legal_d_aggr_thrs

  function is_legal_s_omega(ip)
    implicit none
    real(psb_spk_), intent(in) :: ip
    logical             :: is_legal_s_omega
    is_legal_s_omega = ((ip>=0.0).and.(ip<=2.0))
    return
  end function is_legal_s_omega
  function is_legal_s_fact_thrs(ip)
    implicit none
    real(psb_spk_), intent(in) :: ip
    logical             :: is_legal_s_fact_thrs

    is_legal_s_fact_thrs = (ip>=0.0)
    return
  end function is_legal_s_fact_thrs
  function is_legal_s_aggr_thrs(ip)
    implicit none
    real(psb_spk_), intent(in) :: ip
    logical             :: is_legal_s_aggr_thrs

    is_legal_s_aggr_thrs = (ip>=0.0)
    return
  end function is_legal_s_aggr_thrs


  subroutine amg_icheck_def(ip,name,id,is_legal)
    implicit none
    integer(psb_ipk_), intent(inout) :: ip
    integer(psb_ipk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface
      function is_legal(i)
        import :: psb_ipk_
        integer(psb_ipk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface
    character(len=20), parameter :: rname='amg_check_def'

    if (.not.is_legal(ip)) then
      write(0,*)trim(rname),': Error: Illegal value for ',&
           & name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine amg_icheck_def

  subroutine amg_scheck_def(ip,name,id,is_legal)
    implicit none
    real(psb_spk_), intent(inout) :: ip
    real(psb_spk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface
      function is_legal(i)
        use psb_base_mod, only : psb_spk_
        real(psb_spk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface
    character(len=20), parameter :: rname='amg_check_def'

    if (.not.is_legal(ip)) then
      write(0,*)trim(rname),': Error: Illegal value for ',&
           & name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine amg_scheck_def

  subroutine amg_dcheck_def(ip,name,id,is_legal)
    implicit none
    real(psb_dpk_), intent(inout) :: ip
    real(psb_dpk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface
      function is_legal(i)
        use psb_base_mod, only : psb_dpk_
        real(psb_dpk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface
    character(len=20), parameter :: rname='amg_check_def'

    if (.not.is_legal(ip)) then
      write(0,*)trim(rname),': Error: Illegal value for ',&
           & name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine amg_dcheck_def


  function pr_to_str(iprec)
    implicit none

    integer(psb_ipk_), intent(in)  :: iprec
    character(len=10)     :: pr_to_str

    select case(iprec)
    case(amg_noprec_)
      pr_to_str='NOPREC'
    case(amg_jac_)
      pr_to_str='JAC'
    case(amg_bjac_)
      pr_to_str='BJAC'
    case(amg_as_)
      pr_to_str='AS'
    end select

  end function pr_to_str

  subroutine amg_ml_bcast(ctxt,dat,root)

    implicit none
    type(psb_ctxt_type), intent(in) :: ctxt
    type(amg_ml_parms), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    call psb_bcast(ctxt,dat%sweeps_pre,root)
    call psb_bcast(ctxt,dat%sweeps_post,root)
    call psb_bcast(ctxt,dat%ml_cycle,root)
    call psb_bcast(ctxt,dat%aggr_type,root)
    call psb_bcast(ctxt,dat%par_aggr_alg,root)
    call psb_bcast(ctxt,dat%aggr_ord,root)
    call psb_bcast(ctxt,dat%aggr_prol,root)
    call psb_bcast(ctxt,dat%aggr_omega_alg,root)
    call psb_bcast(ctxt,dat%aggr_eig,root)
    call psb_bcast(ctxt,dat%aggr_filter,root)
    call psb_bcast(ctxt,dat%coarse_mat,root)
    call psb_bcast(ctxt,dat%coarse_solve,root)

  end subroutine amg_ml_bcast

  subroutine amg_sml_bcast(ctxt,dat,root)

    implicit none
    type(psb_ctxt_type), intent(in) :: ctxt
    type(amg_sml_parms), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    call psb_bcast(ctxt,dat%amg_ml_parms,root)
    call psb_bcast(ctxt,dat%aggr_omega_val,root)
    call psb_bcast(ctxt,dat%aggr_thresh,root)
  end subroutine amg_sml_bcast

  subroutine amg_dml_bcast(ctxt,dat,root)
    implicit none
    type(psb_ctxt_type), intent(in) :: ctxt
    type(amg_dml_parms), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    call psb_bcast(ctxt,dat%amg_ml_parms,root)
    call psb_bcast(ctxt,dat%aggr_omega_val,root)
    call psb_bcast(ctxt,dat%aggr_thresh,root)
  end subroutine amg_dml_bcast

  subroutine ml_parms_clone(pm,pmout,info)

    implicit none
    class(amg_ml_parms), intent(inout) :: pm
    class(amg_ml_parms), intent(out)   :: pmout
    integer(psb_ipk_), intent(out)     :: info

    info = psb_success_
    pmout%sweeps_pre     = pm%sweeps_pre
    pmout%sweeps_post    = pm%sweeps_post
    pmout%ml_cycle       = pm%ml_cycle
    pmout%aggr_type      = pm%aggr_type
    pmout%par_aggr_alg   = pm%par_aggr_alg
    pmout%aggr_ord       = pm%aggr_ord
    pmout%aggr_prol      = pm%aggr_prol
    pmout%aggr_omega_alg = pm%aggr_omega_alg
    pmout%aggr_eig       = pm%aggr_eig
    pmout%aggr_filter    = pm%aggr_filter
    pmout%coarse_mat     = pm%coarse_mat
    pmout%coarse_solve   = pm%coarse_solve

  end subroutine ml_parms_clone

  subroutine s_ml_parms_clone(pm,pmout,info)

    implicit none
    class(amg_sml_parms), intent(inout) :: pm
    class(amg_ml_parms), intent(out)   :: pmout
    integer(psb_ipk_), intent(out)     :: info


    integer(psb_ipk_) :: err_act
    integer(psb_ipk_) :: ierr(5)
    character(len=20)  :: name='clone'

    info = 0
    select type(pout => pmout)
    class is (amg_sml_parms)
      call pm%amg_ml_parms%clone(pout%amg_ml_parms,info)
      pout%aggr_omega_val = pm%aggr_omega_val
      pout%aggr_thresh    = pm%aggr_thresh
    class default
      info = psb_err_invalid_dynamic_type_
      ierr(1) = 2
      info = psb_err_missing_override_method_
      call psb_errpush(info,name,i_err=ierr)
      call psb_get_erraction(err_act)
      call psb_error_handler(err_act)
    end select

  end subroutine s_ml_parms_clone

  subroutine d_ml_parms_clone(pm,pmout,info)

    implicit none
    class(amg_dml_parms), intent(inout) :: pm
    class(amg_ml_parms), intent(out)   :: pmout
    integer(psb_ipk_), intent(out)     :: info


    integer(psb_ipk_) :: err_act
    integer(psb_ipk_) :: ierr(5)
    character(len=20)  :: name='clone'

    info = 0
    select type(pout => pmout)
    class is (amg_dml_parms)
      call pm%amg_ml_parms%clone(pout%amg_ml_parms,info)
      pout%aggr_omega_val = pm%aggr_omega_val
      pout%aggr_thresh    = pm%aggr_thresh
    class default
      info = psb_err_invalid_dynamic_type_
      ierr(1) = 2
      info = psb_err_missing_override_method_
      call psb_errpush(info,name,i_err=ierr)
      call psb_get_erraction(err_act)
      call psb_error_handler(err_act)
      return
    end select

  end subroutine d_ml_parms_clone

  function amg_s_equal_aggregation(parms1, parms2) result(val)
    type(amg_sml_parms), intent(in) :: parms1, parms2
    logical :: val

    val  = (parms1%par_aggr_alg     == parms2%par_aggr_alg        ) .and. &
         & (parms1%aggr_type        == parms2%aggr_type       ) .and. &
         & (parms1%aggr_ord         == parms2%aggr_ord        ) .and. &
         & (parms1%aggr_prol        == parms2%aggr_prol       ) .and. &
         & (parms1%aggr_omega_alg   == parms2%aggr_omega_alg  ) .and. &
         & (parms1%aggr_eig         == parms2%aggr_eig        ) .and. &
         & (parms1%aggr_filter      == parms2%aggr_filter     ) .and. &
         & (parms1%aggr_omega_val   == parms2%aggr_omega_val  ) .and. &
         & (parms1%aggr_thresh      == parms2%aggr_thresh     )
  end function amg_s_equal_aggregation

  function amg_d_equal_aggregation(parms1, parms2) result(val)
    type(amg_dml_parms), intent(in) :: parms1, parms2
    logical :: val

    val  = (parms1%par_aggr_alg     == parms2%par_aggr_alg        ) .and. &
         & (parms1%aggr_type        == parms2%aggr_type       ) .and. &
         & (parms1%aggr_ord         == parms2%aggr_ord        ) .and. &
         & (parms1%aggr_prol        == parms2%aggr_prol       ) .and. &
         & (parms1%aggr_omega_alg   == parms2%aggr_omega_alg  ) .and. &
         & (parms1%aggr_eig         == parms2%aggr_eig        ) .and. &
         & (parms1%aggr_filter      == parms2%aggr_filter     ) .and. &
         & (parms1%aggr_omega_val   == parms2%aggr_omega_val  ) .and. &
         & (parms1%aggr_thresh      == parms2%aggr_thresh     )
  end function amg_d_equal_aggregation

  subroutine i_ag_default(ag)
    class(amg_iaggr_data), intent(inout) :: ag

    ag%min_coarse_size = -ione
    ag%min_coarse_size_per_process = -ione
    ag%max_levs    = 20_psb_ipk_
  end subroutine i_ag_default

  subroutine s_ag_default(ag)
    class(amg_saggr_data), intent(inout) :: ag

    call ag%amg_iaggr_data%default()
    ag%min_cr_ratio   = 1.5_psb_spk_
    ag%op_complexity  = szero
    ag%avg_cr         = szero
  end subroutine s_ag_default

  subroutine d_ag_default(ag)
    class(amg_daggr_data), intent(inout) :: ag

    call ag%amg_iaggr_data%default()
    ag%min_cr_ratio   = 1.5_psb_dpk_
    ag%op_complexity  = dzero
    ag%avg_cr         = dzero
  end subroutine d_ag_default


end module amg_base_prec_type

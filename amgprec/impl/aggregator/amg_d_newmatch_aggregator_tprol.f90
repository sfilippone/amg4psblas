!   
!  
! File: amg_d_newmatch_aggregator_tprol.f90
!
! Subroutine: amg_d_newmatch_aggregator_tprol
! Version:    real
!
!
  
subroutine  amg_d_newmatch_aggregator_build_tprol(ag,parms,ag_data,&
     & a,desc_a,ilaggr,nlaggr,op_prol,info)
  use psb_base_mod
  use amg_d_prec_type
  use amg_d_newmatch_aggregator_mod, amg_protect_name => amg_d_newmatch_aggregator_build_tprol
  use amg_d_inner_mod
  use iso_c_binding      
  implicit none
  class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
  type(amg_dml_parms), intent(inout)    :: parms 
  type(amg_daggr_data), intent(in)      :: ag_data
  type(psb_dspmat_type), intent(inout)  :: a
  type(psb_desc_type), intent(inout)    :: desc_a
  integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  type(psb_ldspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)      :: info


  ! Local variables
  real(psb_dpk_), allocatable:: valaggr(:)
  type(psb_dspmat_type)   :: a_tmp
  type(nwm_CSRMatrix) :: C, P
  integer(c_int) :: match_algorithm, n_sweeps, max_csize, max_nlevels
  character(len=20)            :: name, ch_err
  type(psb_ctxt_type)          :: ctxt
  integer(psb_mpk_)            :: np, me
  integer(psb_ipk_)            :: err_act, ierr
  integer(psb_ipk_)            :: debug_level, debug_unit
  integer(psb_ipk_)            :: i, j, k, nr, nc, isz, num_pcols
  type(psb_d_csr_sparse_mat), target :: acsr
  integer(psb_ipk_), allocatable, target ::  csr_ia(:), csr_ja(:), c_ilaggr(:)
  integer(psb_ipk_), allocatable :: aux(:)
  real(psb_dpk_), allocatable, target::  csr_val(:)
  interface
    function bootCMatch(C,match_alg,n_sweeps,max_nlevels,max_csize,w)&
         & bind(c,name='bootCMatch') result(P)
      use iso_c_binding  
      import
      implicit none
      type(nwm_CSRMatrix) :: C, P
      type(nwm_Vector) :: w
      integer(c_int) :: match_alg
      integer(c_int) :: n_sweeps
      integer(c_int) :: max_nlevels
      integer(c_int) :: max_csize
    end function bootCMatch
  end interface

  interface
    function amg_bootCMatch_if(C,match_alg,n_sweeps,max_nlevels,max_csize,&
         & w,isz,ilaggr,valaggr, num_cols) &
         & bind(c,name='amg_bootCMatch_if') result(iret)
      use iso_c_binding  
      import
      implicit none
      type(nwm_CSRMatrix) :: C, P
      type(nwm_Vector) :: w
      integer(c_int), value :: match_alg
      integer(c_int), value :: n_sweeps
      integer(c_int), value :: max_nlevels
      integer(c_int), value :: max_csize
      integer(c_int), value :: isz
      integer(c_int)        :: num_cols
      integer(c_int)        :: ilaggr(*)
      real(c_double)        :: valaggr(*)
      integer(c_int) :: iret
    end function amg_bootCMatch_if
  end interface

  name='amg_d_newmatch_aggregator_tprol'
  ctxt = desc_a%get_context()
  call psb_info(ctxt,me,np)
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_


  call amg_check_def(parms%ml_cycle,'Multilevel cycle',&
       &   amg_mult_ml_,is_legal_ml_cycle)
  call amg_check_def(parms%par_aggr_alg,'Aggregation',&
       &   amg_dec_aggr_,is_legal_decoupled_par_aggr_alg)
  call amg_check_def(parms%aggr_ord,'Ordering',&
       &   amg_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call amg_check_def(parms%aggr_thresh,'Aggr_Thresh',dzero,is_legal_d_aggr_thrs)

  call a%csclip(b=a_tmp, info=info, jmax=a%get_nrows(), imax=a%get_nrows())

  call a_tmp%mv_to(acsr)
  if (ag%do_clean_zeros) call acsr%clean_zeros(info)
  nr = a%get_nrows()
  if (psb_size(ag%w) < nr) call ag%bld_default_w(nr)
  
  !write(*,*) 'Build_tprol:',acsr%get_nrows(),acsr%get_ncols()
  C%num_rows     = acsr%get_nrows()
  C%num_cols     = acsr%get_ncols()
  C%num_nonzeros = acsr%get_nzeros()
  C%owns_data    = 0
  acsr%irp = acsr%irp - 1
  acsr%ja  = acsr%ja  - 1
  C%i    = c_loc(acsr%irp)
  C%j    = c_loc(acsr%ja)
  C%data = c_loc(acsr%val)

  isz = a%get_ncols()
  call psb_realloc(isz,ilaggr,info)
  if (info == psb_success_) call psb_realloc(isz,c_ilaggr,info)
  if (info == psb_success_) call psb_realloc(isz,valaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  match_algorithm = ag%matching_alg
  n_sweeps        = ag%n_sweeps
  if (ag%max_csize > 0) then
    max_csize       = ag%max_csize
  else
    max_csize       = (ag_data%min_coarse_size + np -1)/np
  end if
  if (ag%max_nlevels > 0) then
    max_nlevels     = ag%max_nlevels
  else
    max_nlevels = ag_data%max_levs
  end if

  info = amg_bootCMatch_if(C,match_algorithm,n_sweeps,max_nlevels,max_csize,&
       & ag%w_c_nxt, isz, c_ilaggr, valaggr, num_pcols)
  if (info /= psb_success_) then
!!$      write(0,*) 'On return from bootCMatch_if:',info
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_bootCMatch_if')
    goto 9999
  end if
  ilaggr(1:nr) = c_ilaggr(1:nr)
!!$  write(0,*) 'On output from BootCMatch',nr,num_pcols,size(ilaggr),maxval(ilaggr),&
!!$       & minval(ilaggr),minval(ilaggr(1:nr)),a%get_nrows(),a%get_ncols()
  ! Prepare vector W for next level, just in case
  call ag%bld_wnxt(ilaggr(1:nr),valaggr(1:nr),num_pcols)
  
  call psb_realloc(np,nlaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/np,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if
  call acsr%free()

  nlaggr(:)=0
  nlaggr(me+1) = num_pcols
  call psb_sum(ctxt,nlaggr(1:np))


  call amg_d_newmatch_map_to_tprol(desc_a,ilaggr,nlaggr,valaggr,op_prol,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_newmatch_map_to_tprol')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_newmatch_aggregator_build_tprol

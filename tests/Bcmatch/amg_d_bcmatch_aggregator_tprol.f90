!  
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 , 2017 
!  
!                        Salvatore Filippone  Cranfield University
!  		      Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
! File: mld_d_bcmatch_aggregator_tprol.f90
!
! Subroutine: mld_d_bcmatch_aggregator_tprol
! Version:    real
!
!
  
subroutine  mld_d_bcmatch_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
  use psb_base_mod
  use mld_d_prec_type
  use mld_d_bcmatch_aggregator_mod, mld_protect_name => mld_d_bcmatch_aggregator_build_tprol
  use mld_d_inner_mod
  use bcm_csr_type_mod
  use iso_c_binding      
  implicit none
  class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
  type(mld_dml_parms), intent(inout)  :: parms 
  type(psb_dspmat_type), intent(in)   :: a
  type(psb_desc_type), intent(in)     :: desc_a
  integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  type(psb_dspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)      :: info


  ! Local variables
  real(psb_dpk_), allocatable:: valaggr(:)
  type(psb_dspmat_type)   :: a_tmp
  type(bcm_CSRMatrix) :: C, P
  integer(c_int) :: match_algorithm, n_sweeps, max_csize, max_nlevels
  character(len=20)           :: name, ch_err
  integer(psb_mpk_)           :: ictxt, np, me
  integer(psb_ipk_)           :: err_act, ierr
  integer(psb_ipk_)           :: debug_level, debug_unit
  integer(psb_ipk_)           :: i, j, k, nr, nc, isz, num_pcols
  type(psb_d_csr_sparse_mat), target :: acsr
  integer(psb_ipk_), allocatable, target ::  csr_ia(:), csr_ja(:)
  integer(psb_ipk_), allocatable :: aux(:)
  real(psb_dpk_), allocatable, target::  csr_val(:)
  interface
    function bootCMatch(C,match_alg,n_sweeps,max_nlevels,max_csize,w)bind(c,name='bootCMatch') result(P)
      use iso_c_binding  
      use bcm_csr_type_mod         
      implicit none
      type(bcm_CSRMatrix) :: C, P
      type(bcm_Vector) :: w
      integer(c_int) :: match_alg
      integer(c_int) :: n_sweeps
      integer(c_int) :: max_nlevels
      integer(c_int) :: max_csize
    end function bootCMatch
  end interface

  interface
    function mld_bootCMatch_if(C,match_alg,n_sweeps,max_nlevels,max_csize,&
         & w,isz,ilaggr,valaggr, num_cols) &
         & bind(c,name='mld_bootCMatch_if') result(iret)
      use iso_c_binding  
      use bcm_csr_type_mod         
      implicit none
      type(bcm_CSRMatrix) :: C, P
      type(bcm_Vector) :: w
      integer(c_int), value :: match_alg
      integer(c_int), value :: n_sweeps
      integer(c_int), value :: max_nlevels
      integer(c_int), value :: max_csize
      integer(c_int), value :: isz
      integer(c_int)        :: num_cols
      integer(c_int)        :: ilaggr(*)
      real(c_double)        :: valaggr(*)
      integer(c_int) :: iret
    end function mld_bootCMatch_if
  end interface

  name='mld_d_bcmatch_aggregator_tprol'
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)
  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_


  call mld_check_def(parms%ml_cycle,'Multilevel cycle',&
       &   mld_mult_ml_,is_legal_ml_cycle)
  call mld_check_def(parms%par_aggr_alg,'Aggregation',&
       &   mld_dec_aggr_,is_legal_ml_par_aggr_alg)
  call mld_check_def(parms%aggr_ord,'Ordering',&
       &   mld_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call mld_check_def(parms%aggr_thresh,'Aggr_Thresh',dzero,is_legal_d_aggr_thrs)

  call a%csclip(b=a_tmp, info=info, jmax=a%get_nrows(), imax=a%get_nrows())

  call a_tmp%mv_to(acsr)
  nr = a%get_nrows()
  if (psb_size(ag%w_tmp) < nr) call ag%bld_default_w(nr)
  
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
  if (info == psb_success_) call psb_realloc(isz,valaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  match_algorithm = ag%matching_alg
  n_sweeps        = ag%n_sweeps
  max_csize       = ag%max_csize
  max_nlevels     = ag%max_nlevels

  info = mld_bootCMatch_if(C,match_algorithm,n_sweeps,max_nlevels,max_csize,&
       & ag%w_par, isz, ilaggr, valaggr, num_pcols)
  if (info /= psb_success_) then
!!$      write(0,*) 'On return from bootCMatch_if:',info
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_bootCMatch_if')
    goto 9999
  end if
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
  call psb_sum(ictxt,nlaggr(1:np))


  call mld_d_bcmatch_map_to_tprol(desc_a,ilaggr,nlaggr,valaggr,op_prol,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mld_bcmatch_map_to_tprol')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine mld_d_bcmatch_aggregator_build_tprol

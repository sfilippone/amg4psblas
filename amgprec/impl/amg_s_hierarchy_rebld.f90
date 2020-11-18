!   
!   
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
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
! File: amg_s_hierarchy_bld.f90
!
! Subroutine: amg_s_hierarchy_bld
! Version:    real
!
!  This routine builds the preconditioner according to the requirements made by
!  the user trough the subroutines amg_precinit and amg_precset.
!  
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level,
!  (for more details see the description of amg_Tonelev_type in amg_prec_type.f90).
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level. No transfer operators are associated to level 1.
! 
!
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_dprec_type), input/output.
!               The preconditioner data structure; upon exit it contains 
!               the multilevel hierarchy of prolongators, restrictors
!               and coarse matrices.
!    info    -  integer, output.
!               Error code.              
!  
subroutine amg_s_hierarchy_rebld(a,desc_a,prec,info)

  use psb_base_mod
  use amg_s_inner_mod
  use amg_s_prec_mod, amg_protect_name => amg_s_hierarchy_rebld

  Implicit None

  ! Arguments
  type(psb_sspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  class(amg_sprec_type),intent(inout),target          :: prec
  integer(psb_ipk_), intent(out)                       :: info

  ! Local Variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  integer(psb_ipk_)   :: err,i,k, err_act, iszv, newsz,&
       & nplevs, mxplevs
  integer(psb_lpk_) :: iaggsize, casize
  real(psb_spk_)     :: mnaggratio, sizeratio, athresh, aomega
  class(amg_s_base_smoother_type), allocatable :: coarse_sm, med_sm, &
       & med_sm2, coarse_sm2
  class(amg_s_base_aggregator_type), allocatable :: tmp_aggr
  type(amg_sml_parms)              :: medparms, coarseparms
  integer(psb_lpk_), allocatable   :: nlaggr(:)
  type(psb_s_coo_sparse_mat) :: coo_prol, coo_restr
  type(psb_s_csr_sparse_mat) :: acsr
  type(psb_desc_type), pointer :: p_desc_a
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'amg_hierarchy_rebld'
  info = psb_success_
  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  prec%ctxt = ctxt
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  !

  if (.not.allocated(prec%precv)) then 
    !! Error: should have called amg_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if


  iszv = size(prec%precv)

  do i=2, iszv
    call prec%precv(i-1)%base_a%cp_to(acsr)
    p_desc_a => prec%precv(i-1)%base_desc
    call prec%precv(i)%linmap%mat_V2U%cp_to(coo_prol)
    call prec%precv(i)%linmap%mat_U2V%cp_to(coo_restr)
    call amg_rap(acsr,p_desc_a,prec%precv(i)%linmap%naggr,&
         & prec%precv(i)%parms,prec%precv(i)%ac,&
         & coo_prol,prec%precv(i)%desc_ac,coo_restr,info)

  end do

  call prec%cmp_complexity()
  call prec%cmp_avg_cr()

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Exiting with',iszv,' levels'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_s_hierarchy_rebld

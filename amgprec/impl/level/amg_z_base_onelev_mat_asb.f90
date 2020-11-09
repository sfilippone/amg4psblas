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
! File: amg_z_onelev_mat_asb.f90
!
! Subroutine: amg_z_onelev_mat_asb
! Version:    complex
!
!  This routine builds the matrix associated to the current level of the
!  multilevel preconditioner from the matrix associated to the previous level,
!  by using the user-specified aggregation technique (therefore, it also builds the
!  prolongation and restriction operators mapping the current level to the
!  previous one and vice versa). 
!  The current level is regarded as the coarse one, while the previous as
!  the fine one. This is in agreement with the fact that the routine is called,
!  by amg_mlprec_bld, only on levels >=2.
!  The main structure is:
!  1. Perform sanity checks;
!  2. Call amg_Xaggrmat_asb to compute prolongator/restrictor/AC
!  3. According to the choice of DIST/REPL for AC, build a descriptor DESC_AC,
!     and adjust the column numbering of AC/OP_PROL/OP_RESTR
!  4. Pack restrictor and prolongator into p%linmap
!  5. Fix base_a and base_desc pointers.
!
! 
! Arguments:
!    p       -  type(amg_z_onelev_type), input/output.
!               The 'one-level' data structure containing the control
!               parameters and (eventually) coarse matrix and prolongator/restrictors. 
!               
!    a       -  type(psb_zspmat_type).
!               The sparse matrix structure containing the local part of the
!               fine-level matrix.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    ilaggr     -  integer, dimension(:), input
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix. Note that the indices
!                  are assumed to be shifted so as to make sure the ranges on
!                  the various processes do not   overlap.
!    nlaggr     -  integer, dimension(:) input
!                  nlaggr(i) contains the aggregates held by process i.
!    op_prol    -  type(psb_zspmat_type), input/output
!               The tentative prolongator on input, released on output. 
!               
!    info    -  integer, output.
!               Error code.         
!  
subroutine amg_z_base_onelev_mat_asb(lv,a,desc_a,ilaggr,nlaggr,t_prol,info)

  use psb_base_mod
  use amg_base_prec_type
  use amg_z_onelev_mod, amg_protect_name => amg_z_base_onelev_mat_asb

  implicit none

  ! Arguments
  class(amg_z_onelev_type), intent(inout), target :: lv
  type(psb_zspmat_type), intent(in)  :: a
  type(psb_desc_type), intent(inout)   :: desc_a
  integer(psb_lpk_), intent(inout) :: nlaggr(:)
  integer(psb_lpk_), intent(inout) :: ilaggr(:)
  type(psb_lzspmat_type), intent(inout)  :: t_prol
  integer(psb_ipk_), intent(out)      :: info
  

  ! Local variables
  character(len=24)            :: name
  integer(psb_ipk_)            :: ictxt, np, me
  integer(psb_ipk_)            :: err_act
  type(psb_zspmat_type)        :: ac, op_restr, op_prol
  integer(psb_ipk_)            :: nzl, inl
  integer(psb_ipk_)            :: debug_level, debug_unit

  name='amg_z_onelev_mat_asb'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  call amg_check_def(lv%parms%aggr_prol,'Smoother',&
       &   amg_smooth_prol_,is_legal_ml_aggr_prol)
  call amg_check_def(lv%parms%coarse_mat,'Coarse matrix',&
       &   amg_distr_mat_,is_legal_ml_coarse_mat)
  call amg_check_def(lv%parms%aggr_filter,'Use filtered matrix',&
       &   amg_no_filter_mat_,is_legal_aggr_filter)
  call amg_check_def(lv%parms%aggr_omega_alg,'Omega Alg.',&
       &   amg_eig_est_,is_legal_ml_aggr_omega_alg)
  call amg_check_def(lv%parms%aggr_eig,'Eigenvalue estimate',&
       &   amg_max_norm_,is_legal_ml_aggr_eig)
  call amg_check_def(lv%parms%aggr_omega_val,'Omega',dzero,is_legal_d_omega)


  !
  ! Build the coarse-level matrix from the fine-level one, starting from 
  ! the mapping defined by amg_aggrmap_bld and applying the aggregation
  ! algorithm specified by lv%iprcparm(amg_aggr_prol_)
  !
  call lv%aggr%mat_bld(lv%parms,a,desc_a,ilaggr,nlaggr,&
       & lv%ac,lv%desc_ac,op_prol,op_restr,t_prol,info)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_aggrmat_asb')
    goto 9999
  end if

  !
  ! Now build its descriptor and convert global indices for
  ! ac, op_restr and op_prol
  !
  if (info == psb_success_) &
       & call lv%aggr%mat_asb(lv%parms,a,desc_a,&
       & lv%ac,lv%desc_ac,op_prol,op_restr,info)
  
  if (info == psb_success_) call lv%ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
  
  if (info == psb_success_) call lv%aggr%bld_map(desc_a, lv%desc_ac,&
       & ilaggr,nlaggr,op_restr,op_prol,lv%linmap,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mat_asb/map_bld')
    goto 9999
  end if
  !
  ! Fix the base_a and base_desc pointers for handling of residuals.
  ! This is correct because this routine is only called at levels >=2.
  !
  lv%base_a    => lv%ac
  lv%base_desc => lv%desc_ac

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_z_base_onelev_mat_asb

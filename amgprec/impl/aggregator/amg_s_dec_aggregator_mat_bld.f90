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
! File: amg_s_dec_aggregator_mat_bld.f90
!
! Subroutine: amg_s_dec_aggregator_mat_bld
! Version:    real
!
!  This routine builds the matrix associated to the current level of the
!  multilevel preconditioner from the matrix associated to the previous level,
!  by using the user-specified aggregation technique (therefore, it also builds the
!  prolongation and restriction operators mapping the current level to the
!  previous one and vice versa). 
!  The current level is regarded as the coarse one, while the previous as
!  the fine one. This is in agreement with the fact that the routine is called,
!  by amg_mlprec_bld, only on levels >=2.
!  The coarse-level matrix A_C is built from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
! 
!  A mapping from the nodes of the adjacency graph of A to the nodes of the
!  adjacency graph of A_C has been computed by the amg_aggrmap_bld subroutine.
!  The prolongator P_C is built here from this mapping, according to the
!  value of p%iprcparm(amg_aggr_kind_), specified by the user through
!  amg_sprecinit and amg_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  amg_s_lev_aggrmat_bld.
!
!  Currently four  different prolongators are implemented, corresponding to
!  four  aggregation algorithms:
!  1. un-smoothed aggregation,
!  2. smoothed aggregation,
!  3. "bizarre" aggregation.
!  4. minimum energy 
!  1. The non-smoothed aggregation uses as prolongator the piecewise constant
!     interpolation operator corresponding to the fine-to-coarse level mapping built
!     by p%aggr%bld_tprol. This is called tentative prolongator.
!  2. The smoothed aggregation uses as prolongator the operator obtained by applying
!     a damped Jacobi smoother to the tentative prolongator.
!  3. The "bizarre" aggregation uses a prolongator proposed by the authors of MLD2P4.
!     This prolongator still requires a deep analysis and testing and its use is
!     not recommended.
!  4. Minimum energy aggregation
!
!  For more details see
!    M. Brezina and P. Vanek, A black-box iterative solver based on a two-level
!    Schwarz method, Computing,  63 (1999), 233-263.
!    P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
!    (1996), 179-196.
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of PSBLAS-based
!    parallel two-level Schwarz preconditioners, Appl. Num. Math., 57 (2007),
!    1181-1196.
!    M. Sala, R. Tuminaro: A new Petrov-Galerkin smoothed aggregation preconditioner
!    for nonsymmetric linear systems, SIAM J. Sci. Comput., 31(1):143-166 (2008)
!
!
!  The main structure is:
!  1. Perform sanity checks;
!  2. Compute prolongator/restrictor/AC
!
! 
! Arguments:
!    ag       -  type(amg_s_dec_aggregator_type), input/output.
!               The aggregator object
!    parms   -  type(amg_sml_parms), input 
!               The aggregation parameters
!    a          -  type(psb_sspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
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
!    ac         -  type(psb_sspmat_type), output
!                  The coarse matrix on output 
!                  
!    op_prol    -  type(psb_sspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_sspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!  
subroutine  amg_s_dec_aggregator_mat_bld(ag,parms,a,desc_a,ilaggr,nlaggr,&
     & ac,desc_ac,op_prol,op_restr,t_prol,info)
  use psb_base_mod
  use amg_s_prec_type, amg_protect_name => amg_s_dec_aggregator_mat_bld
  use amg_s_inner_mod 
  implicit none
  
  class(amg_s_dec_aggregator_type), target, intent(inout) :: ag
  type(amg_sml_parms), intent(inout)      :: parms 
  type(psb_sspmat_type), intent(in)    :: a
  type(psb_desc_type), intent(inout)     :: desc_a
  integer(psb_lpk_), intent(inout)     :: ilaggr(:), nlaggr(:)
  type(psb_lsspmat_type), intent(inout)   :: t_prol
  type(psb_sspmat_type), intent(out)   :: op_prol,ac,op_restr
  type(psb_desc_type), intent(inout)       :: desc_ac
  integer(psb_ipk_), intent(out)       :: info

  ! Local variables
  character(len=20)            :: name
  integer(psb_mpk_)            :: ictxt, np, me
  integer(psb_lpk_)            :: nzl,ntaggr
  integer(psb_ipk_)            :: err_act
  integer(psb_ipk_)            :: debug_level, debug_unit

  name='amg_s_dec_aggregator_mat_bld'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  !
  ! Build the coarse-level matrix from the fine-level one, starting from 
  ! the mapping defined by amg_aggrmap_bld and applying the aggregation
  ! algorithm specified by 
  !
  select case (parms%aggr_prol)
  case (amg_no_smooth_) 

    call amg_saggrmat_nosmth_bld(a,desc_a,ilaggr,nlaggr,&
         & parms,ac,desc_ac,op_prol,op_restr,t_prol,info)

  case(amg_smooth_prol_) 

    call amg_saggrmat_smth_bld(a,desc_a,ilaggr,nlaggr, &
         & parms,ac,desc_ac,op_prol,op_restr,t_prol,info)

!!$  case(amg_biz_prol_) 
!!$
!!$    call amg_saggrmat_biz_bld(a,desc_a,ilaggr,nlaggr, &
!!$         & parms,ac,desc_ac,op_prol,op_restr,t_prol,info)

  case(amg_min_energy_) 

    call amg_saggrmat_minnrg_bld(a,desc_a,ilaggr,nlaggr, &
         & parms,ac,desc_ac,op_prol,op_restr,t_prol,info)

  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid aggr kind')
    goto 9999

  end select
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner aggrmat bld')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

  
end subroutine amg_s_dec_aggregator_mat_bld

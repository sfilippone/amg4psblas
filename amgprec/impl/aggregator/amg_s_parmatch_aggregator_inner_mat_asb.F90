!
!
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.7)
!
!    (C) Copyright 2021
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
! File: amg_s_parmatch_aggregator_mat_asb.f90
!
! Subroutine: amg_s_parmatch_aggregator_mat_asb
! Version:    real
!
!
!  From a given AC to final format, generating DESC_AC.
!  This is quite involved, because in the context of aggregation based
!  on parallel matching we are building the matrix hierarchy within BLD_TPROL
!  as we go, especially if we have multiple sweeps, hence this code is called
!  in two completely different contexts:
!  1. Within bld_tprol for the internal hierarchy
!  2. Outside, from amg_hierarchy_bld
!  The solution we have found is for bld_tprol to copy its output
!  into special components ag%ac ag%desc_ac etc so that:
!  1. if they are  allocated, it means that bld_tprol has been already invoked, we are in
!     amg_hierarchy_bld and we only need to copy them
!  2. If they are not allocated, we are within bld_tprol, and we need to actually
!     perform the various needed steps.
!
! Arguments:
!    ag       -  type(amg_s_parmatch_aggregator_type), input/output.
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
!    ac         -  type(psb_sspmat_type), inout
!                  The coarse matrix
!    desc_ac    -  type(psb_desc_type), output.
!                  The communication descriptor of the fine-level matrix.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!
!    op_prol    -  type(psb_sspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!
!    op_restr    -  type(psb_sspmat_type), input/output
!                  The restrictor operator; normally, it is the transpose of the prolongator.
!
!    info       -  integer, output.
!                  Error code.
!
subroutine  amg_s_parmatch_aggregator_inner_mat_asb(ag,parms,a,desc_a,&
     & ac,desc_ac, op_prol,op_restr,info)
  use psb_base_mod
  use amg_base_prec_type
#if defined(SERIAL_MPI)
    use amg_s_parmatch_aggregator_mod
#else
  use amg_s_parmatch_aggregator_mod, amg_protect_name => amg_s_parmatch_aggregator_inner_mat_asb
#endif
  implicit none
  class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
  type(amg_sml_parms), intent(inout)    :: parms
  type(psb_sspmat_type), intent(in)     :: a
  type(psb_desc_type), intent(in)       :: desc_a
  type(psb_sspmat_type), intent(inout) :: op_prol,op_restr
  type(psb_sspmat_type), intent(inout)  :: ac
  type(psb_desc_type), intent(inout)    :: desc_ac
  integer(psb_ipk_), intent(out)        :: info
  !
  type(psb_ctxt_type)         :: ictxt
  integer(psb_ipk_)           :: np, me
  type(psb_ls_coo_sparse_mat) :: acoo, bcoo
  type(psb_ls_csr_sparse_mat) :: acsr1
  integer(psb_ipk_)           :: nzl, inl
  integer(psb_lpk_)           :: ntaggr
  integer(psb_ipk_) :: err_act, debug_level, debug_unit
  character(len=20) :: name='d_parmatch_inner_mat_asb'
  character(len=80) :: aname
  logical, parameter :: debug=.false., dump_prol_restr=.false.


  if (psb_get_errstatus().ne.0) return
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

#if !defined(SERIAL_MPI)

  if (debug) write(0,*) me,' ',trim(name),' Start:',&
       & allocated(ag%ac),allocated(ag%desc_ac), allocated(ag%prol),allocated(ag%restr)

  select case(parms%coarse_mat)

  case(amg_distr_mat_)
    ! Do nothing, it has already been done in spmm_bld_ov.

  case(amg_repl_mat_)
    !
    !
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='no repl coarse_mat_ here')
    goto 9999

  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid amg_coarse_mat_')
    goto 9999
  end select
#endif
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return


end subroutine amg_s_parmatch_aggregator_inner_mat_asb

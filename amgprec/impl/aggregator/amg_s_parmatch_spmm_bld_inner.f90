!
!
!                             AMG4PSBLAS  Extensions
!
!    (C) Copyright 2019
!
!                        Salvatore Filippone  Cranfield University
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
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
! File: amg_saggrmat_nosmth_bld.F90
!
! Subroutine: amg_saggrmat_nosmth_bld
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is the piecewise constant interpolation operator corresponding
!  the fine-to-coarse level mapping built by amg_aggrmap_bld.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat
!  specified by the user through amg_sprecinit and amg_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    a          -  type(psb_sspmat_type), input.
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(amg_s_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(amg_sml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_sspmat_type), output
!                  The coarse matrix on output
!
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
!    op_prol    -  type(psb_sspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!
!    op_restr    -  type(psb_sspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator.
!
!    info       -  integer, output.
!                  Error code.
!
!
subroutine amg_s_parmatch_spmm_bld_inner(a_csr,desc_a,ilaggr,nlaggr,parms,&
     & ac,desc_ac,op_prol,op_restr,t_prol,info)
  use psb_base_mod
  use amg_s_inner_mod
  use amg_s_parmatch_aggregator_mod, amg_protect_name => amg_s_parmatch_spmm_bld_inner
  implicit none

  ! Arguments
  type(psb_s_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_lpk_), intent(inout)          :: ilaggr(:), nlaggr(:)
  type(amg_sml_parms), intent(inout)        :: parms
  type(psb_lsspmat_type), intent(inout)     :: t_prol
  type(psb_sspmat_type), intent(out)        :: ac, op_prol, op_restr
  type(psb_desc_type), intent(out)          :: desc_ac
  integer(psb_ipk_), intent(out)            :: info

  ! Local variables
  integer(psb_ipk_)   :: err_act
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_)   :: np, me, ndx
  character(len=40)   :: name
  type(psb_ls_coo_sparse_mat) :: tmpcoo
  type(psb_s_coo_sparse_mat)  :: coo_prol, coo_restr
  type(psb_s_csr_sparse_mat)  :: ac_csr, csr_restr
  type(psb_desc_type), target  :: tmp_desc
  type(psb_lsspmat_type)     :: lac
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nrl, nzl, ip, &
       &  nzt, naggrm1, naggrp1, i, k
  integer(psb_lpk_), allocatable :: ia(:),ja(:)
  !integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza, nrpsave, ncpsave, nzpsave
  logical, parameter :: do_timings=.true., oldstyle=.false., debug=.false.
  integer(psb_ipk_), save :: idx_spspmm=-1, idx_prolcnv=-1, idx_proltrans=-1, idx_asb=-1

  name='amg_parmatch_spmm_bld_inner'
  if(psb_get_errstatus().ne.0) return
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("SPMM_BLD: spspmm   ")
  if ((do_timings).and.(idx_prolcnv==-1)) &
       & idx_prolcnv = psb_get_timer_idx("SPMM_BLD: prolcnv  ")
  if ((do_timings).and.(idx_proltrans==-1)) &
       & idx_proltrans = psb_get_timer_idx("SPMM_BLD: proltrans")
  if ((do_timings).and.(idx_asb==-1)) &
       & idx_asb = psb_get_timer_idx("SPMM_BLD: asb      ")

  if (do_timings) call psb_tic(idx_prolcnv)
  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))

  !
  ! Here T_PROL should be arriving with GLOBAL indices on the cols
  ! and LOCAL indices on the rows.
  !
  if (debug) write(0,*)  me,' ',trim(name),' Size check on entry New: ',&
       & op_prol%get_fmt(),op_prol%get_nrows(),op_prol%get_ncols(),op_prol%get_nzeros(),&
       & nrow,ntaggr,naggr

  call t_prol%cp_to(tmpcoo)

  call psb_cdall(ictxt,desc_ac,info,nl=naggr)
  nzl = tmpcoo%get_nzeros()
  if (debug) write(0,*)  me,' ',trim(name),' coo_prol: ',&
       & tmpcoo%ia(1:min(10,nzl)),' :',tmpcoo%ja(1:min(10,nzl))
  call desc_ac%indxmap%g2lip_ins(tmpcoo%ja(1:nzl),info)
  call tmpcoo%set_ncols(desc_ac%get_local_cols())
  call tmpcoo%cp_to_icoo(coo_prol,info)

  call amg_ptap_bld(a_csr,desc_a,nlaggr,parms,ac,&
       & coo_prol,desc_ac,coo_restr,info)

  nzl = coo_prol%get_nzeros()
  if (debug) write(0,*)  me,' ',trim(name),' coo_prol: ',&
       & coo_prol%ia(1:min(10,nzl)),' :',coo_prol%ja(1:min(10,nzl))

  call op_prol%mv_from(coo_prol)
  call op_restr%mv_from(coo_restr)

  if (debug) then
    write(0,*) me,' ',trim(name),' Checkpoint at exit'
    call psb_barrier(ictxt)
    write(0,*) me,' ',trim(name),' Checkpoint through'
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build ac = op_restr x a3')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done smooth_aggregate '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_s_parmatch_spmm_bld_inner

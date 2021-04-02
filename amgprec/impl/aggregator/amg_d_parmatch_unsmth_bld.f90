!
!
!                             AMG4PSBLAS  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!
!    (C) Copyright 2008-2018
!
!        Salvatore Filippone
!        Pasqua D'Ambra
!        Daniela di Serafino
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
! File: amg_d_parmatch_unsmth_bld.F90
!
! Subroutine: amg_d_parmatch_unsmth_bld
! Version:    real
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
!
!  The prolongator P_C is built according to a smoothed aggregation algorithm,
!  i.e. it is obtained by applying a damped Jacobi smoother to the piecewise
!  constant interpolation operator P corresponding to the fine-to-coarse level
!  mapping built by the amg_aggrmap_bld subroutine:
!
!                            P_C = (I - omega*D^(-1)A) * P,
!
!  where D is the diagonal matrix with main diagonal equal to the main diagonal
!  of A, and omega is a suitable smoothing parameter. An estimate of the spectral
!  radius of D^(-1)A, to be used in the computation of omega, is provided,
!  according to the value of p%parms%aggr_omega_alg, specified by the user
!  through amg_dprecinit and amg_zprecset.
!
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat,
!  specified by the user through amg_dprecinit and amg_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  aggregator%mat_bld.
!
!
! Arguments:
!    a          -  type(psb_dspmat_type), input.
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(amg_d_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(amg_dml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_dspmat_type), output
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
!    op_prol    -  type(psb_dspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!
!    op_restr    -  type(psb_dspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator.
!
!    info       -  integer, output.
!                  Error code.
!
subroutine amg_d_parmatch_unsmth_bld(ag,a,desc_a,ilaggr,nlaggr,parms,&
     & ac,desc_ac,op_prol,op_restr,t_prol,info)
  use psb_base_mod
  use amg_base_prec_type
  use amg_d_inner_mod
  use amg_d_base_aggregator_mod
  use amg_d_parmatch_aggregator_mod, amg_protect_name => amg_d_parmatch_unsmth_bld
  implicit none

  ! Arguments
  class(amg_d_parmatch_aggregator_type), target, intent(inout) :: ag
  type(psb_dspmat_type), intent(in)      :: a
  type(psb_desc_type), intent(inout)     :: desc_a
  integer(psb_lpk_), intent(inout)         :: ilaggr(:), nlaggr(:)
  type(amg_dml_parms), intent(inout)    :: parms
  type(psb_dspmat_type), intent(inout)    :: op_prol,ac,op_restr
  type(psb_ldspmat_type), intent(inout)  :: t_prol
  type(psb_desc_type), intent(inout)       :: desc_ac
  integer(psb_ipk_), intent(out)           :: info

  ! Local variables
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, ip, &
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrw
  integer(psb_ipk_) :: inaggr
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_)   :: np, me
  character(len=20)   :: name
  type(psb_ld_coo_sparse_mat) :: lcoo_prol
  type(psb_d_coo_sparse_mat) :: coo_prol, coo_restr
  type(psb_d_csr_sparse_mat) :: acsr
  type(psb_d_csr_sparse_mat)  :: csr_prol, acsr3, csr_restr, ac_csr
  real(psb_dpk_), allocatable :: adiag(:)
  real(psb_dpk_), allocatable :: arwsum(:)
  logical            :: filter_mat
  integer(psb_ipk_)            :: debug_level, debug_unit, err_act
  integer(psb_ipk_), parameter :: ncmax=16
  real(psb_dpk_)   :: anorm, omega, tmp, dg, theta
  logical, parameter :: debug_new=.false., dump_r=.false., dump_p=.false., debug=.false.
  logical, parameter :: do_timings=.false.
  integer(psb_ipk_), save :: idx_spspmm=-1
  character(len=80) :: filename

  name='amg_parmatch_unsmth_bld'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  theta = parms%aggr_thresh
  !write(0,*) me,' ',trim(name),' Start '

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("PMC_UNSMTH_BLD: par_spspmm")

  !
  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))
  !write(0,*) me,' ',trim(name),' input sizes',nlaggr(:),':',naggr

  call a%cp_to(acsr)
  call t_prol%mv_to(lcoo_prol)

  inaggr = naggr
  call psb_cdall(ictxt,desc_ac,info,nl=inaggr)
  nzl = lcoo_prol%get_nzeros()
  call desc_ac%indxmap%g2lip_ins(lcoo_prol%ja(1:nzl),info)
  call lcoo_prol%set_ncols(desc_ac%get_local_cols())
  call lcoo_prol%cp_to_icoo(coo_prol,info)

  if (debug) call check_coo(me,trim(name)//' Check 1 on  coo_prol:',coo_prol)

  call psb_cdasb(desc_ac,info)
  call psb_cd_reinit(desc_ac,info)
  if (.not.allocated(ag%desc_ax)) allocate(ag%desc_ax)

  call amg_ptap_bld(acsr,desc_a,nlaggr,parms,ac,&
       & coo_prol,desc_ac,coo_restr,info,desc_ax=ag%desc_ax)

  call op_restr%cp_from(coo_restr)
  call op_prol%mv_from(coo_prol)

  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr
  ! write(0,*)  me,' ',trim(name),' Final AC newstyle ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros()

  if (debug) then
    write(0,*) me,' ',trim(name),' Checkpoint at exit'
    call psb_barrier(ictxt)
    write(0,*) me,' ',trim(name),' Checkpoint through'
    block
      character(len=128) :: fname, prefix_
      integer :: lname
      prefix_ = "unsmth_bld_"
      lname = len_trim(prefix_)
      fname = trim(prefix_)
      write(fname(lname+1:lname+10),'(a,i3.3,a)') '_p_',me, '.mtx'
      call op_prol%print(fname,head='Debug aggregates')
      write(fname(lname+1:lname+10),'(a,i3.3,a)') '_r_',me, '.mtx'
      call op_restr%print(fname,head='Debug aggregates')
      write(fname(lname+1:lname+11),'(a,i3.3,a)') '_ac_',me, '.mtx'
      call ac%print(fname,head='Debug aggregates')
    end block
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Build ac = coo_restr x am3')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return

contains
  subroutine check_coo(me,string,coo)
    implicit none
    integer(psb_ipk_) :: me
    type(psb_d_coo_sparse_mat) :: coo
    character(len=*) :: string
    integer(psb_lpk_) :: nr,nc,nz
    nr = coo%get_nrows()
    nc = coo%get_ncols()
    nz = coo%get_nzeros()
    write(0,*) me,string,nr,nc,&
         & minval(coo%ia(1:nz)),maxval(coo%ia(1:nz)),&
         & minval(coo%ja(1:nz)),maxval(coo%ja(1:nz))
  end subroutine check_coo

end subroutine amg_d_parmatch_unsmth_bld

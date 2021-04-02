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
! File: amg_daggrmat_smth_bld.F90
!
! Subroutine: amg_daggrmat_smth_bld
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
subroutine amg_d_parmatch_smth_bld(ag,a,desc_a,ilaggr,nlaggr,parms,&
     & ac,desc_ac,op_prol,op_restr,t_prol,info)
  use psb_base_mod
  use amg_base_prec_type
  use amg_d_inner_mod
  use amg_d_base_aggregator_mod
  use amg_d_parmatch_aggregator_mod, amg_protect_name => amg_d_parmatch_smth_bld
  implicit none

  ! Arguments
  class(amg_d_parmatch_aggregator_type), target, intent(inout) :: ag
  type(psb_dspmat_type), intent(in)      :: a
  type(psb_desc_type), intent(inout)     :: desc_a
  integer(psb_lpk_), intent(inout)         :: ilaggr(:), nlaggr(:)
  type(amg_dml_parms), intent(inout)    :: parms
  type(psb_ldspmat_type), intent(inout)  :: t_prol
  type(psb_dspmat_type), intent(out)    :: op_prol,ac,op_restr
  type(psb_desc_type), intent(inout)    :: desc_ac
  integer(psb_ipk_), intent(out)           :: info

  ! Local variables
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, ip, &
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrw
  integer(psb_ipk_)   :: inaggr
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_)   :: np, me
  character(len=20)   :: name
  type(psb_ld_coo_sparse_mat) :: tmpcoo, ac_coo, lcoo_restr
  type(psb_d_coo_sparse_mat)  :: coo_prol, coo_restr
  type(psb_d_csr_sparse_mat)  :: acsrf, csr_prol, acsr, tcsr
  real(psb_dpk_), allocatable :: adiag(:)
  real(psb_dpk_), allocatable :: arwsum(:)
  logical            :: filter_mat
  integer(psb_ipk_)            :: debug_level, debug_unit, err_act
  integer(psb_ipk_), parameter :: ncmax=16
  real(psb_dpk_)     :: anorm, omega, tmp, dg, theta
  logical, parameter :: debug_new=.false., dump_r=.false., dump_p=.false., debug=.false.
  character(len=80) :: filename
  logical, parameter :: do_timings=.false.
  integer(psb_ipk_), save :: idx_spspmm=-1, idx_phase1=-1, idx_gtrans=-1, idx_phase2=-1, idx_refine=-1, idx_phase3=-1
  integer(psb_ipk_), save :: idx_cdasb=-1, idx_ptap=-1

  name='amg_parmatch_smth_bld'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !debug_level = 2
  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  theta = parms%aggr_thresh
  !write(0,*) me,' ',trim(name),' Start ',idx_spspmm
  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("PMC_SMTH_BLD: par_spspmm")
  if ((do_timings).and.(idx_phase1==-1)) &
       & idx_phase1 = psb_get_timer_idx("PMC_SMTH_BLD: phase1    ")
  if ((do_timings).and.(idx_phase2==-1)) &
       & idx_phase2 = psb_get_timer_idx("PMC_SMTH_BLD: phase2    ")
  if ((do_timings).and.(idx_phase3==-1)) &
       & idx_phase3 = psb_get_timer_idx("PMC_SMTH_BLD: phase3    ")
  if ((do_timings).and.(idx_gtrans==-1)) &
       & idx_gtrans = psb_get_timer_idx("PMC_SMTH_BLD: gtrans    ")
  if ((do_timings).and.(idx_refine==-1)) &
       & idx_refine = psb_get_timer_idx("PMC_SMTH_BLD: refine    ")
  if ((do_timings).and.(idx_cdasb==-1)) &
       & idx_cdasb = psb_get_timer_idx("PMC_SMTH_BLD: cdasb     ")
  if ((do_timings).and.(idx_ptap==-1)) &
       & idx_ptap = psb_get_timer_idx("PMC_SMTH_BLD: ptap_bld  ")

  if (do_timings) call psb_tic(idx_phase1)

  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))
  filter_mat = (parms%aggr_filter == amg_filter_mat_)

  !
  ! naggr: number of local aggregates
  ! nrow: local rows.
  !
  if (dump_p) then
    block
      integer(psb_lpk_), allocatable :: ivr(:), ivc(:)
      integer(psb_lpk_) :: i
      character(len=132) :: aname
      type(psb_ldspmat_type) :: aglob
      type(psb_dspmat_type) :: atmp
!!$      call a%cp_to(acsr)
!!$      call atmp%cp_from(acsr)
      write(0,*) me,' ',trim(name),' Dumping inp_prol/restr'
      write(aname,'(a,i0,a,i0,a)') 'tprol-',desc_a%get_global_rows(),'-p',me,'.mtx'
      call t_prol%print(fname=aname,head='Test ')
    end block
  end if

  if (do_timings)  call psb_tic(idx_refine)
  ! Get the diagonal D
  adiag = a%get_diag(info)
  if (info == psb_success_) &
       & call psb_realloc(ncol,adiag,info)
  if (info == psb_success_) &
       & call psb_halo(adiag,desc_a,info)
  if (info == psb_success_) call a%cp_to(acsr)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Initial copies done.'

  call acsr%cp_to_fmt(acsrf,info)

  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    !

    do i=1, nrow
      tmp = dzero
      jd  = -1
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) jd = j
        if (abs(acsrf%val(j)) < theta*sqrt(abs(adiag(i)*adiag(acsrf%ja(j))))) then
          tmp=tmp+acsrf%val(j)
          acsrf%val(j)=dzero
        endif

      enddo
      if (jd == -1) then
        write(0,*) 'Wrong input: we need the diagonal!!!!', i
      else
        acsrf%val(jd)=acsrf%val(jd)-tmp
      end if
    enddo
    ! Take out zeroed terms
    call acsrf%clean_zeros(info)
  end if


  do i=1,size(adiag)
    if (adiag(i) /= dzero) then
      adiag(i) = done / adiag(i)
    else
      adiag(i) = done
    end if
  end do
  if (do_timings) call psb_toc(idx_refine)

  if (parms%aggr_omega_alg == amg_eig_est_) then

    if (parms%aggr_eig == amg_max_norm_) then
      allocate(arwsum(nrow))
      call acsr%arwsum(arwsum)
      anorm = maxval(abs(adiag(1:nrow)*arwsum(1:nrow)))
      call psb_amx(ictxt,anorm)
      omega = 4.d0/(3.d0*anorm)
      parms%aggr_omega_val = omega

    else
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid amg_aggr_eig_')
      goto 9999
    end if

  else if (parms%aggr_omega_alg == amg_user_choice_) then

    omega = parms%aggr_omega_val

  else if (parms%aggr_omega_alg /= amg_user_choice_) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid amg_aggr_omega_alg_')
    goto 9999
  end if


  call acsrf%scal(adiag,info)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Filtering and scaling  done.',info
  if (info /= psb_success_) goto 9999

  inaggr = naggr

  call t_prol%cp_to(tmpcoo)

  call psb_cdall(ictxt,desc_ac,info,nl=inaggr)
  nzl = tmpcoo%get_nzeros()
  call desc_ac%indxmap%g2lip_ins(tmpcoo%ja(1:nzl),info)
  call tmpcoo%set_ncols(desc_ac%get_local_cols())
  call tmpcoo%mv_to_ifmt(tcsr,info)
  !
  ! Build the smoothed prolongator using either A or Af
  !    csr_prol = (I-w*D*A) Prol      csr_prol = (I-w*D*Af) Prol
  ! This is always done through the variable acsrf which
  ! is a bit less readable, but saves space and one extra matrix copy
  !
  call omega_smooth(omega,acsrf)
  if (do_timings) call psb_toc(idx_phase1)

  if (do_timings) call psb_tic(idx_spspmm)
  call psb_par_spspmm(acsrf,desc_a,tcsr,csr_prol,desc_ac,info)
  call tcsr%free()
  if (do_timings) call psb_toc(idx_spspmm)
  if (do_timings) call psb_tic(idx_phase2)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 1')
    goto 9999
  end if
  !
  ! Now that we have the smoothed prolongator, we can
  ! compute the triple product.
  !
  if (do_timings)  call psb_tic(idx_cdasb)
  call psb_cdasb(desc_ac,info)
  if (do_timings)  call psb_toc(idx_cdasb)
  call psb_cd_reinit(desc_ac,info)

  call csr_prol%mv_to_coo(coo_prol,info)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done SPSPMM 1'

  if (do_timings)  call psb_tic(idx_ptap)
  if (.not.allocated(ag%desc_ax)) allocate(ag%desc_ax)
  call amg_ptap_bld(acsr,desc_a,nlaggr,parms,ac,&
       & coo_prol,desc_ac,coo_restr,info,desc_ax=ag%desc_ax)
  if (do_timings)  call psb_toc(idx_ptap)

  call op_prol%mv_from(coo_prol)
  call op_restr%mv_from(coo_restr)


  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr
  ! write(0,*)  me,' ',trim(name),' Final AC newstyle ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros()

  if (dump_r) then
    block
      integer(psb_lpk_), allocatable :: ivr(:), ivc(:)
      integer(psb_lpk_) :: i
      character(len=132) :: aname
      type(psb_ldspmat_type) :: aglob
      type(psb_dspmat_type) :: atmp
      write(0,*) me,' ',trim(name),' Dumping prol/restr'
      ivc    = [(i,i=1,desc_a%get_local_cols())]
      call desc_a%l2gip(ivc,info)
      ivr    = [(i,i=1,desc_ac%get_local_cols())]
      call desc_ac%l2gip(ivr,info)

      write(aname,'(a,i0,a,i0,a)') 'restr-',desc_ac%get_global_rows(),'-p',me,'.mtx'

      call op_restr%print(fname=aname,head='Test ',ivc=ivc)
!!$      write(aname,'(a,i0,a,i0,a)') 'prol-',desc_ac%get_global_rows(),'-p',me,'.mtx'
!!$      call op_prol%print(fname=aname,head='Test ')
!!$      call psb_gather(aglob,atmp,desc_a,info)
!!$      if (me==psb_root_) then
!!$        write(aname,'(a,i0,a)') 'a-inp-g-',aglob%get_nrows(),'.mtx'
!!$        call aglob%print(fname=aname,head='Test ')
!!$      end if

    end block
  end if
  if (do_timings) call psb_toc(idx_phase2)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done smooth_aggregate '
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return

contains

  subroutine omega_smooth(omega,acsr)
    implicit none
    real(psb_dpk_),intent(in) :: omega
    type(psb_d_csr_sparse_mat), intent(inout) :: acsr
    !
    integer(psb_ipk_) :: i,j
    do i=1,acsr%get_nrows()
      do j=acsr%irp(i),acsr%irp(i+1)-1
        if (acsr%ja(j) == i) then
          acsr%val(j) = done - omega*acsr%val(j)
        else
          acsr%val(j) = - omega*acsr%val(j)
        end if
      end do
    end do
  end subroutine omega_smooth

end subroutine amg_d_parmatch_smth_bld

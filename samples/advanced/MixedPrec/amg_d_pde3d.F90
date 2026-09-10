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
!         software without specific prior written permission.
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
!
! File: amg_d_pde3d.f90
!
! Program: amg_d_pde3d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs.
!
!
! The PDE is a general second order equation in 3d
!
!   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)
! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
!      dxdx     dydy       dzdz        dx       dy         dz
!
! with Dirichlet boundary conditions
!   u = g
!
!  on the unit cube  0<=x,y,z<=1.
!
!
! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
!
! There are three choices available for data distribution:
! 1. A simple BLOCK distribution
! 2. A ditribution based on arbitrary assignment of indices to processes,
!    typically from a graph partitioner
! 3. A 3D distribution in which the unit cube is partitioned
!    into subcubes, each one assigned to a process.
!
!
! File: amg_d_pde3d.F90
!
! Module: psb_cg_mixed_mod
!
!   This module contains a copy of the PSBLAS psb_dcg implementation,
!   renamed psb_cg_mixed and specialized for mixed precision: the Krylov
!   recurrence, the matrix-vector products and the convergence tests are
!   carried out in double precision, whereas the preconditioner is built
!   and applied in single precision.
!
!   At each iteration the application of the preconditioner
!
!               z = M^(-1) r
!
!   is replaced by
!
!               r_single = single(r)
!               z_single = M^(-1) r_single      (single precision preconditioner)
!               z  = double(z_single)
!
!   the conversions being performed by the psb_d2s_vect/psb_s2d_vect
!   routines of psb_mixed_support_mod. Everything else is identical to
!   psb_dcg.
!
module psb_cg_mixed_mod
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod
  use psb_mixed_support_mod

  implicit none

contains

  !
  ! Subroutine: psb_cg_mixed
  !    Conjugate Gradient in double precision with a single precision
  !    preconditioner.
  !
  ! Arguments:
  !
  !    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A,
  !                                         in double precision.
  !    prec   -  class(psb_sprec_type)      Input: preconditioner, in SINGLE
  !                                         precision; it must have been built
  !                                         on the single precision copy of A.
  !    b      -  type(psb_d_vect_type)      Input: right hand side B.
  !    x      -  type(psb_d_vect_type)      Input/Output: initial guess and
  !                                         final solution X.
  !    eps    -  real                       Input: stopping tolerance; the
  !                                         iteration is stopped when the error
  !                                         estimate |err| <= eps
  !    desc_a -  type(psb_desc_type)        Input: the communication descriptor;
  !                                         it is precision independent, hence
  !                                         it is shared by the double and the
  !                                         single precision vectors.
  !    info   -  integer                    Output: return code
  !
  !    itmax  -  integer(optional)          Input: maximum number of iterations
  !                                         to be performed.
  !    iter   -  integer(optional)          Output: how many iterations have
  !                                         been performed.
  !    err    -  real   (optional)          Output: error estimate on exit.
  !    itrace -  integer(optional)          Input: print an informational message
  !                                         with the error estimate every itrace
  !                                         iterations
  !    istop  -  integer(optional)          Input: stopping criterion.
  !    cond   -  real   (optional)          Output: condition number estimate.
  !
  subroutine psb_cg_mixed(a,prec,b,x,eps,desc_a,info,&
       & itmax,iter,err,itrace,istop,cond)
    implicit none
    type(psb_dspmat_type), intent(in)    :: a
    Type(psb_desc_type), Intent(in)      :: desc_a
    class(psb_sprec_type), intent(inout) :: prec
    type(psb_d_vect_type), Intent(inout) :: b
    type(psb_d_vect_type), Intent(inout) :: x
    Real(psb_dpk_), Intent(in)           :: eps
    integer(psb_ipk_), intent(out)                 :: info
    integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
    integer(psb_ipk_), Optional, Intent(out)       :: iter
    Real(psb_dpk_), Optional, Intent(out) :: err,cond
! =   Local data
    real(psb_dpk_), allocatable, target   :: aux(:),td(:),tu(:),eig(:),ewrk(:)
    integer(psb_mpk_), allocatable :: ibl(:), ispl(:), iwrk(:)
    type(psb_d_vect_type), allocatable, target :: wwrk(:)
    type(psb_d_vect_type), pointer  :: q, p, r, z, w
! =   Single precision workspace for the preconditioner application
    real(psb_spk_), allocatable, target   :: saux(:)
    type(psb_s_vect_type)                 :: r_single, z_single
    real(psb_dpk_)   :: alpha, beta, rho, rho_old, sigma,alpha_old,beta_old
    integer(psb_ipk_) :: itmax_, istop_, naux, it, itx, itrace_,&
         &  n_col, n_row,err_act, ieg,nspl, istebz
    integer(psb_lpk_) :: mglob
    integer(psb_ipk_) :: debug_level, debug_unit
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: np, me
    real(psb_dpk_)     :: derr
    type(psb_itconv_type)       :: stopdat
    logical                     :: do_cond
    character(len=20)           :: name
    character(len=*), parameter :: methdname='CG-MIXED'

    info = psb_success_
    name = 'psb_cg_mixed'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ctxt = desc_a%get_context()

    call psb_info(ctxt, me, np)
    if (.not.allocated(b%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
    endif


    mglob = desc_a%get_global_rows()
    n_row = desc_a%get_local_rows()
    n_col = desc_a%get_local_cols()


    if (present(istop)) then
      istop_ = istop
    else
      istop_ = psb_get_istop_default()
    endif
    if (.not.psb_is_valid_istop(istop_)) then
      info=psb_err_invalid_istop_
      err=info
      call psb_errpush(info,name,i_err=(/istop_/))
      goto 9999
    end if
    !
    !  istop_ = 1:  normwise backward error, infinity norm
    !  istop_ = 2:  ||r||/||b||   norm 2
    !
    select case(istop_)
    case(psb_istop_ani_,psb_istop_bn2_,&
         & psb_istop_rn2_abs_, psb_istop_rrn2_)
      ! nothing needed
    case default
      ! should never get here
      info=psb_err_internal_error_
      err=info
      call psb_errpush(info,name,a_err="invalid istop_")
      goto 9999
    end select

    call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
    if (info == psb_success_)&
         & call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on X/B')
      goto 9999
    end if

    naux=4*n_col
    allocate(aux(naux),saux(naux), stat=info)
    if (info == psb_success_) call psb_geall(wwrk,desc_a,info,n=5_psb_ipk_)
    if (info == psb_success_) call psb_geasb(wwrk,desc_a,info,mold=x%v,scratch=.true.)
    !
    !  The descriptor is precision independent: the single precision
    !  vectors handed to the preconditioner are allocated on the very
    !  same descriptor as the double precision ones.
    !
    if (info == psb_success_) call psb_geall(r_single,desc_a,info)
    if (info == psb_success_) call psb_geasb(r_single,desc_a,info,scratch=.true.)
    if (info == psb_success_) call psb_geall(z_single,desc_a,info)
    if (info == psb_success_) call psb_geasb(z_single,desc_a,info,scratch=.true.)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

    p  => wwrk(1)
    q  => wwrk(2)
    r  => wwrk(3)
    z  => wwrk(4)
    w  => wwrk(5)
    !
    !  The vectors crossing the precision boundary are converted in
    !  full, halo included; the halo of R is only written by the
    !  preconditioner, hence the scratch storage is cleaned here to
    !  avoid converting uninitialized values.
    !
    call r%zero()
    call z_single%zero()


    if (present(itmax)) then
      itmax_ = itmax
    else
      itmax_ = 1000
    endif

    if (present(itrace)) then
      itrace_ = itrace
    else
      itrace_ = 0
    end if

    do_cond=present(cond)
    if (do_cond) then
      istebz = 0
      allocate(td(itmax_),tu(itmax_), eig(itmax_),&
           & ibl(itmax_),ispl(itmax_),iwrk(3*itmax_),ewrk(4*itmax_),&
           & stat=info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
    end if
    itx=0
    alpha = dzero

    restart: do
! =
! =    r0 = b-Ax0
! =
      if (itx>= itmax_) exit restart

      it = 0
      call psb_geaxpby(done,b,dzero,r,desc_a,info)
      if (info == psb_success_) call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if

      rho = dzero

      call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
      if (info /= psb_success_) Then
        call psb_errpush(psb_err_from_subroutine_non_,name)
        goto 9999
      End If

      iteration:  do

        it   = it + 1
        itx = itx + 1

        !
        !  This is the only step performed in single precision:
        !  convert the residual, apply the single precision
        !  preconditioner, then convert the result back to double
        !  and carry on with the double precision recurrence.
        !
        call psb_d2s_vect(r,r_single,info)
        if (info == psb_success_) call prec%apply(r_single,z_single,desc_a,info,work=saux)
        if (info == psb_success_) call psb_s2d_vect(z_single,z,info,mold=x%v)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_non_
          call psb_errpush(info,name,a_err='single precision preconditioner')
          goto 9999
        end if

        rho_old = rho
        rho     = psb_gedot(r,z,desc_a,info)

        if (it == 1) then
          call psb_geaxpby(done,z,dzero,p,desc_a,info)
        else
          if (rho_old == dzero) then
            if (debug_level >= psb_debug_ext_)&
                 & write(debug_unit,*) me,' ',trim(name),&
                 & ': CG Iteration breakdown rho'
            exit iteration
          endif
          beta = rho/rho_old
          call psb_geaxpby(done,z,beta,p,desc_a,info)
        end if

        call psb_spmm(done,a,p,dzero,q,desc_a,info,work=aux)
        sigma = psb_gedot(p,q,desc_a,info)
        if (sigma == dzero) then
          if (debug_level >= psb_debug_ext_)&
               & write(debug_unit,*) me,' ',trim(name),&
               & ': CG Iteration breakdown sigma'
          exit iteration
        endif
        alpha_old = alpha
        alpha = rho/sigma

        if (do_cond) then
          istebz = istebz + 1
          if (istebz == 1) then
            td(istebz) = done/alpha
          else
            td(istebz) = done/alpha + beta/alpha_old
            tu(istebz-1) = sqrt(beta)/alpha_old
          end if
        end if


        call psb_geaxpby(alpha,p,done,x,desc_a,info)
        call psb_geaxpby(-alpha,q,done,r,desc_a,info)

        if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
        if (info /= psb_success_) Then
          call psb_errpush(psb_err_from_subroutine_non_,name)
          goto 9999
        End If

      end do iteration
    end do restart
    if (do_cond) then
      if (me == psb_root_) then
#if defined(PSB_HAVE_LAPACK)
        call dstebz('A','E',istebz,dzero,dzero,0,0,-done,td,tu,&
             & ieg,nspl,eig,ibl,ispl,ewrk,iwrk,info)
        if (info < 0) then
          call psb_errpush(psb_err_from_subroutine_ai_,name,&
               & a_err='dstebz',i_err=(/info/))
          info=psb_err_from_subroutine_ai_
          goto 9999
        end if
        cond = eig(ieg)/eig(1)
#else
        cond = dzero
#endif
        info=psb_success_
      end if
      call psb_bcast(ctxt,cond)
    end if


    call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
    if (present(err)) err = derr

    if (info == psb_success_) call psb_gefree(wwrk,desc_a,info)
    if (info == psb_success_) call psb_gefree(r_single,desc_a,info)
    if (info == psb_success_) call psb_gefree(z_single,desc_a,info)
    if (info == psb_success_) deallocate(aux,saux,stat=info)
    if (info /= psb_success_) then
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_cg_mixed

end module psb_cg_mixed_mod

program amg_d_pde3d
  use psb_base_mod
  use amg_prec_mod
  use psb_linsolve_mod
  use psb_util_mod
  use data_input
  use amg_d_pde3d_poisson_mod
  use amg_d_pde3d_exp_mod
  use amg_d_pde3d_box_mod
  use amg_d_pde3d_gauss_mod
  use amg_d_genpde_mod
  use psb_mixed_support_mod
  use psb_cg_mixed_mod
#if defined(PSB_OPENMP)
  use omp_lib
#endif
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  character(len=32)  :: pdecoeff
  integer(psb_ipk_) :: idim
  integer(psb_epk_) :: system_size

  ! miscellaneous
  real(psb_dpk_) :: t1, t2, tprec, thier, tslv, tsmth, tpgen

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(psb_sspmat_type) :: asingle
  type(amg_dprec_type)  :: prec
  type(amg_sprec_type)  :: sprec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: x,b,r
  ! parallel environment
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np, nth

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize, vecsize
  real(psb_dpk_)   :: err, resmx, resmxp

  ! Solver data
  type solverdata
    character(len=40)  :: kmethd      ! Iterative solver
    integer(psb_ipk_)  :: istopc      ! stopping criterion
    integer(psb_ipk_)  :: itmax       ! maximum number of iterations
    integer(psb_ipk_)  :: itrace      ! tracing
    integer(psb_ipk_)  :: irst        ! restart
    real(psb_dpk_)     :: eps         ! stopping tolerance
  end type solverdata
  type(solverdata)       :: s_choice

  ! preconditioner data
  type precdata

    ! preconditioner type
    character(len=40)  :: descr       ! verbose description of the prec
    character(len=10)  :: ptype       ! preconditioner type

    integer(psb_ipk_)  :: outer_sweeps ! number of outer sweeps: sweeps for 1-level,
    ! AMG cycles for ML
    ! general AMG data
    character(len=32)  :: mlcycle      ! AMG cycle type
    integer(psb_ipk_)  :: maxlevs      ! maximum number of levels in AMG preconditioner

    ! AMG aggregation
    character(len=32)  :: aggr_prol    ! aggregation type: SMOOTHED, NONSMOOTHED
    character(len=32)  :: par_aggr_alg ! parallel aggregation algorithm: DEC, SYMDEC
    character(len=32)  :: aggr_type    ! Type of aggregation SOC1, SOC2, MATCHBOXP
    integer(psb_ipk_)  :: aggr_size    ! Requested size of the aggregates for MATCHBOXP
    character(len=32)  :: aggr_ord     ! ordering for aggregation: NATURAL, DEGREE
    character(len=32)  :: aggr_filter  ! filtering: FILTER, NO_FILTER
    real(psb_dpk_)      :: mncrratio    ! minimum aggregation ratio
    real(psb_dpk_), allocatable :: athresv(:) ! smoothed aggregation threshold vector
    integer(psb_ipk_)  :: thrvsz       ! size of threshold vector
    real(psb_dpk_)      :: athres       ! smoothed aggregation threshold
    integer(psb_ipk_)  :: csizepp      ! minimum size of coarsest matrix per process

    ! AMG smoother or pre-smoother; also 1-lev preconditioner
    character(len=32)  :: smther       ! (pre-)smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps      ! (pre-)smoother / 1-lev prec. sweeps
    integer(psb_ipk_)  :: degree       ! degree for polynomial smoother
    character(len=32)  :: pvariant     ! polynomial  variant
    character(len=32)  :: prhovariant  ! how to estimate rho(M^{-1}A)
    real(psb_dpk_)     :: prhovalue    ! if previous is set value, we set it from this one
    integer(psb_ipk_)  :: novr         ! number of overlap layers
    character(len=32)  :: restr        ! restriction over application of AS
    character(len=32)  :: prol         ! prolongation over application of AS
    character(len=32)  :: solve        ! local subsolver type: ILU, MILU, ILUT,
                                       ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: ssweeps      ! inner solver sweeps
    character(len=32)  :: variant      ! AINV variant: LLK, etc
    integer(psb_ipk_)  :: fill         ! fill-in for incomplete LU factorization
    integer(psb_ipk_)  :: invfill      ! Inverse fill-in for INVK
    real(psb_dpk_)      :: thr          ! threshold for ILUT factorization

    ! AMG post-smoother; ignored by 1-lev preconditioner
    character(len=32)  :: smther2      ! post-smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps2     ! post-smoother sweeps
    integer(psb_ipk_)  :: degree2      ! degree for polynomial smoother
    character(len=32)  :: pvariant2    ! polynomial  variant
    character(len=32)  :: prhovariant2 ! how to estimate rho(M^{-1}A)
    real(psb_dpk_)     :: prhovalue2   ! if previous is set value, we set it from this one
    integer(psb_ipk_)  :: novr2        ! number of overlap layers
    character(len=32)  :: restr2       ! restriction  over application of AS
    character(len=32)  :: prol2        ! prolongation over application of AS
    character(len=32)  :: solve2       ! local subsolver type: ILU, MILU, ILUT,
                                       ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: ssweeps2     ! inner solver sweeps
    character(len=32)  :: variant2     ! AINV variant: LLK, etc
    integer(psb_ipk_)  :: fill2        ! fill-in for incomplete LU factorization
    integer(psb_ipk_)  :: invfill2     ! Inverse fill-in for INVK
    real(psb_dpk_)      :: thr2         ! threshold for ILUT factorization

    ! coarsest-level solver
    character(len=32)  :: cmat         ! coarsest matrix layout: REPL, DIST
    character(len=32)  :: csolve       ! coarsest-lev solver: BJAC, SLUDIST (distr.
                                       ! mat.); UMF, MUMPS, SLU, ILU, ILUT, MILU
                                       ! (repl. mat.)
    character(len=32)  :: csbsolve     ! coarsest-lev local subsolver: ILU, ILUT,
                                       ! MILU, UMF, MUMPS, SLU
    integer(psb_ipk_)  :: cfill        ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: cthres        ! threshold for ILUT factorization
    integer(psb_ipk_)  :: cjswp        ! sweeps for GS or JAC coarsest-lev subsolver
    ! settings for the Krylov method
    character(len=16)  :: krm_method  ! Krylov method for coarsest level
    character(len=16)  :: krm_prec    ! Preconditioner for coarsest level
    character(len=16)  :: krm_subsolve ! Subsolver for coarsest level
    character(len=16)  :: krm_global  ! Is the solver global or local? TRUE or FALSE
    real(psb_dpk_)      :: krm_eps     ! Stopping tolerance
    integer(psb_ipk_)  :: krm_irst    ! Restart for Krylov method 
    integer(psb_ipk_)  :: krm_istop   ! Stopping criterion
    integer(psb_ipk_)  :: krm_itmax   ! Maximum number of iterations
    integer(psb_ipk_)  :: krm_itrace  ! Trace of the lower Krylov iterations
    integer(psb_ipk_)  :: krm_fillin  ! Fill-in for incomplete LU factorization

    ! Dump data
    logical            :: dump = .false.
    integer(psb_ipk_)  :: dlmin        ! Minimum level to dump 
    integer(psb_ipk_)  :: dlmax        ! Maximum level to dump
    logical            :: dump_ac         = .false.
    logical            :: dump_rp         = .false.
    logical            :: dump_tprol      = .false.
    logical            :: dump_smoother   = .false.
    logical            :: dump_solver     = .false.
    logical            :: dump_global_num = .false.

  end type precdata
  type(precdata)       :: p_choice

  ! other variables
  integer(psb_ipk_)  :: info, i, k
  character(len=20)  :: name,ch_err
  ! .true. when the Krylov method uses the single precision preconditioner
  logical            :: mixed

  info=psb_success_


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
#if defined(PSB_OPENMP)
  !$OMP parallel shared(nth)
  !$OMP master
  nth = omp_get_num_threads()
  !$OMP end master
  !$OMP end parallel
#else
  nth = 1
#endif

  if (iam < 0) then
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='amg_d_pde3d'
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then
    write(*,*) 'Welcome to AMG4PSBLAS version: ',amg_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  !
  !  get parameters
  !
  call get_parms(ctxt,afmt,idim,s_choice,p_choice,pdecoeff)

  !
  !  Only the preconditioner actually used by the chosen Krylov method is
  !  built: the mixed method applies the single precision one, all the
  !  other methods the double precision one.
  !
  select case(psb_toupper(trim(s_choice%kmethd)))
  case('CG-MIXED','CGMIXED')
    mixed = .true.
  case default
    mixed = .false.
  end select

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess
  !

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  select case(psb_toupper(trim(pdecoeff)))
  case("POISSON")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
         & a1_poisson,a2_poisson,a3_poisson,&
         & b1_poisson,b2_poisson,b3_poisson,c_poisson,g_poisson,info)
  case("EXP")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
         & a1_exp,a2_exp,a3_exp,&
         & b1_exp,b2_exp,b3_exp,c_exp,g_exp,info)
  case("BOX")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
         & a1_box,a2_box,a3_box,&
         & b1_box,b2_box,b3_box,c_box,g_box,info)
  case("GAUSS")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
         & a1_gauss,a2_gauss,a3_gauss,&
         & b1_gauss,b2_gauss,b3_gauss,c_gauss,g_gauss,info)
  case default
    info=psb_err_from_subroutine_
    ch_err='amg_gen_pdecoeff'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select


  call psb_barrier(ctxt)
  tpgen = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='amg_gen_pde3d'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("PDE Coefficients             : ",a)')pdecoeff
  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')tpgen
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')
  
  if (mixed) then
  !
  !  Convert the matrix to single precision and build the single
  !  precision preconditioner on it.
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime() 
  call psb_d2s_cscnv(a,asingle,info,type=afmt)
  call psb_barrier(ctxt)
  t2 = psb_wtime() -t1
  if (iam == psb_root_) &
       & write(psb_out_unit,'("Conversion D to S       time : ",es12.5)')t2
  


  !
  ! initialize the preconditioner
  !
  call sprec%init(ctxt,p_choice%ptype,info)
  select case(trim(psb_toupper(p_choice%ptype)))
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','L1-JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call sprec%set('smoother_sweeps', p_choice%jsweeps, info)

  case ('BJAC','POLY')
    call sprec%set('smoother_sweeps', p_choice%jsweeps, info)
    call sprec%set('sub_solve',       p_choice%solve,   info)
    call sprec%set('solver_sweeps',   p_choice%ssweeps,   info)
    call sprec%set('poly_degree',     p_choice%degree,    info)
    call sprec%set('poly_variant',    p_choice%pvariant,  info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call sprec%set('mumps_loc_glob','local_solver',info)
    call sprec%set('sub_fillin',      p_choice%fill,    info)
    call sprec%set('sub_iluthrs',     real(p_choice%thr,kind=psb_spk_),     info)

  case('AS')
    call sprec%set('smoother_sweeps', p_choice%jsweeps, info)
    call sprec%set('sub_ovr',         p_choice%novr,    info)
    call sprec%set('sub_restr',       p_choice%restr,   info)
    call sprec%set('sub_prol',        p_choice%prol,    info)
    call sprec%set('sub_solve',       p_choice%solve,   info)
    call sprec%set('solver_sweeps',   p_choice%ssweeps,   info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call sprec%set('mumps_loc_glob','local_solver',info)
    call sprec%set('sub_fillin',      p_choice%fill,    info)
    call sprec%set('sub_iluthrs',     real(p_choice%thr,kind=psb_spk_),     info)

  case ('ML')
    ! multilevel preconditioner

    call sprec%set('ml_cycle',        p_choice%mlcycle,    info)
    call sprec%set('outer_sweeps',    p_choice%outer_sweeps,info)
    if (p_choice%csizepp>0)&
         & call sprec%set('min_coarse_size_per_process', p_choice%csizepp,    info)
    if (p_choice%mncrratio>1)&
         & call sprec%set('min_cr_ratio',   real(p_choice%mncrratio,kind=psb_spk_), info)
    if (p_choice%maxlevs>0)&
         & call sprec%set('max_levs',    p_choice%maxlevs,    info)
    if (p_choice%athres >= dzero) &
         & call sprec%set('aggr_thresh',     real(p_choice%athres,kind=psb_spk_),  info)
    if (p_choice%thrvsz>0) then
      do k=1,min(p_choice%thrvsz,size(sprec%precv)-1)
        call sprec%set('aggr_thresh',     real(p_choice%athresv(k),kind=psb_spk_),  info,ilev=(k+1))
      end do
    end if

    call sprec%set('aggr_prol',       p_choice%aggr_prol,   info)
    call sprec%set('par_aggr_alg',    p_choice%par_aggr_alg,   info)
    call sprec%set('aggr_type',       p_choice%aggr_type, info)
    call sprec%set('aggr_size',       p_choice%aggr_size, info)

    call sprec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call sprec%set('aggr_filter',     p_choice%aggr_filter,info)


    call sprec%set('smoother_type',   p_choice%smther,     info)
    call sprec%set('smoother_sweeps', p_choice%jsweeps,    info)
    call sprec%set('poly_degree',     p_choice%degree,    info)
    call sprec%set('poly_variant',    p_choice%pvariant,  info)
    if (p_choice%prhovalue > dzero ) then
            call sprec%set('poly_rho_ba', real(p_choice%prhovalue,kind=psb_spk_), info)
    else
            call sprec%set('poly_rho_estimate', p_choice%prhovariant, info)
    end if

    select case (psb_toupper(p_choice%smther))
    case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI','L1-FBGS')
      ! do nothing
    case default
      call sprec%set('sub_ovr',         p_choice%novr,       info)
      call sprec%set('sub_restr',       p_choice%restr,      info)
      call sprec%set('sub_prol',        p_choice%prol,       info)
      select case(trim(psb_toupper(p_choice%solve)))
      case('INVK')
        call sprec%set('sub_solve',       p_choice%solve,   info)
      case('INVT')
        call sprec%set('sub_solve',       p_choice%solve,   info)
      case('AINV')
        call sprec%set('sub_solve',       p_choice%solve,   info)
        call sprec%set('ainv_alg', p_choice%variant,   info)
      case default
        call sprec%set('sub_solve',       p_choice%solve,   info)
        if (psb_toupper(p_choice%solve)=='MUMPS') &
             & call sprec%set('mumps_loc_glob','local_solver',info)
      end select
      call sprec%set('solver_sweeps',   p_choice%ssweeps,   info)
      call sprec%set('sub_fillin',      p_choice%fill,       info)
      call sprec%set('inv_fillin',      p_choice%invfill,    info)
      call sprec%set('sub_iluthrs',     real(p_choice%thr,kind=psb_spk_),        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call sprec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call sprec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
      call sprec%set('poly_degree',     p_choice%degree2,   info,pos='post')
      call sprec%set('poly_variant',    p_choice%pvariant2, info,pos='post')
      if (p_choice%prhovalue > dzero ) then
         call sprec%set('poly_rho_ba', real(p_choice%prhovalue2,kind=psb_spk_), info,pos='post')
      else
         call sprec%set('poly_rho_estimate', p_choice%prhovariant2, info,pos='post')
      end if

      select case (psb_toupper(p_choice%smther2))
      case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI','L1-FBGS')
        ! do nothing
      case default
        call sprec%set('sub_ovr',         p_choice%novr2,     info,pos='post')
        call sprec%set('sub_restr',       p_choice%restr2,    info,pos='post')
        call sprec%set('sub_prol',        p_choice%prol2,     info,pos='post')
        select case(trim(psb_toupper(p_choice%solve2)))
        case('INVK')
          call sprec%set('sub_solve',       p_choice%solve2,   info)
        case('INVT')
          call sprec%set('sub_solve',       p_choice%solve2,   info)
        case('AINV')
          call sprec%set('sub_solve',       p_choice%solve2,   info)
          call sprec%set('ainv_alg', p_choice%variant2,   info)
        case default
          call sprec%set('sub_solve',       p_choice%solve2,   info, pos='post')
          if (psb_toupper(p_choice%solve2)=='MUMPS') &
               & call sprec%set('mumps_loc_glob','local_solver',info)
        end select
        call sprec%set('solver_sweeps',   p_choice%ssweeps2,   info,pos='post')
        call sprec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call sprec%set('inv_fillin',      p_choice%invfill2,  info,pos='post')
        call sprec%set('sub_iluthrs',     real(p_choice%thr2,kind=psb_spk_),      info,pos='post')
      end select
    end if

    call sprec%set('coarse_solve',    p_choice%csolve,    info)
    call sprec%set('coarse_mat',      p_choice%cmat,      info)
     ! Set for the case of a KRM solver
    if (psb_toupper(p_choice%csolve) == 'KRM') then
      call sprec%set('krm_method', p_choice%krm_method, info)
      call sprec%set('krm_kprec', p_choice%krm_prec,   info)
      call sprec%set('krm_sub_solve', p_choice%krm_subsolve, info)
      call sprec%set('krm_global', p_choice%krm_global, info)
      call sprec%set('krm_eps',   real(p_choice%krm_eps,kind=psb_spk_),   info)
      call sprec%set('krm_irst',   p_choice%krm_irst,   info)
      call sprec%set('krm_istopc',  p_choice%krm_istop,  info)
      call sprec%set('krm_itmax',  p_choice%krm_itmax,  info)
      call sprec%set('krm_itrace', p_choice%krm_itrace, info)
      call sprec%set('krm_fillin', p_choice%krm_fillin, info)
    else
      if (psb_toupper(p_choice%csolve) == 'BJAC') &
        &  call sprec%set('coarse_subsolve', p_choice%csbsolve,  info)
      call sprec%set('coarse_fillin',   p_choice%cfill,     info)
      call sprec%set('coarse_iluthrs',  real(p_choice%cthres,kind=psb_spk_),    info)
      call sprec%set('coarse_sweeps',   p_choice%cjswp,     info)
    end if
  end select

  ! build the preconditioner
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call sprec%hierarchy_build(asingle,desc_a,info)
  thier = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_hierarchy_bld')
    goto 9999
  end if
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call sprec%smoothers_build(asingle,desc_a,info)
  tsmth = psb_wtime()-t1 
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ctxt, thier)
  call psb_amx(ctxt, tsmth)
  
  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Single Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Single Preconditioner time: ",es12.5)') thier+tsmth,thier,tsmth
    write(psb_out_unit,'(" ")')
  end if

  else

  !
  ! initialize the preconditioner
  !
  call prec%init(ctxt,p_choice%ptype,info)
  select case(trim(psb_toupper(p_choice%ptype)))
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','L1-JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)

  case ('BJAC','POLY')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('solver_sweeps',   p_choice%ssweeps,   info)
    call prec%set('poly_degree',     p_choice%degree,    info)
    call prec%set('poly_variant',    p_choice%pvariant,  info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call prec%set('mumps_loc_glob','local_solver',info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case('AS')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_ovr',         p_choice%novr,    info)
    call prec%set('sub_restr',       p_choice%restr,   info)
    call prec%set('sub_prol',        p_choice%prol,    info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('solver_sweeps',   p_choice%ssweeps,   info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call prec%set('mumps_loc_glob','local_solver',info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case ('ML')
    ! multilevel preconditioner

    call prec%set('ml_cycle',        p_choice%mlcycle,    info)
    call prec%set('outer_sweeps',    p_choice%outer_sweeps,info)
    if (p_choice%csizepp>0)&
         & call prec%set('min_coarse_size_per_process', p_choice%csizepp,    info)
    if (p_choice%mncrratio>1)&
         & call prec%set('min_cr_ratio',   p_choice%mncrratio, info)
    if (p_choice%maxlevs>0)&
         & call prec%set('max_levs',    p_choice%maxlevs,    info)
    if (p_choice%athres >= dzero) &
         & call prec%set('aggr_thresh',     p_choice%athres,  info)
    if (p_choice%thrvsz>0) then
      do k=1,min(p_choice%thrvsz,size(prec%precv)-1)
        call prec%set('aggr_thresh',     p_choice%athresv(k),  info,ilev=(k+1))
      end do
    end if

    call prec%set('aggr_prol',       p_choice%aggr_prol,   info)
    call prec%set('par_aggr_alg',    p_choice%par_aggr_alg,   info)
    call prec%set('aggr_type',       p_choice%aggr_type, info)
    call prec%set('aggr_size',       p_choice%aggr_size, info)

    call prec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call prec%set('aggr_filter',     p_choice%aggr_filter,info)


    call prec%set('smoother_type',   p_choice%smther,     info)
    call prec%set('smoother_sweeps', p_choice%jsweeps,    info)
    call prec%set('poly_degree',     p_choice%degree,    info)
    call prec%set('poly_variant',    p_choice%pvariant,  info)
    if (p_choice%prhovalue > dzero ) then
            call prec%set('poly_rho_ba', p_choice%prhovalue, info)
    else
            call prec%set('poly_rho_estimate', p_choice%prhovariant, info)
    end if

    select case (psb_toupper(p_choice%smther))
    case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI','L1-FBGS')
      ! do nothing
    case default
      call prec%set('sub_ovr',         p_choice%novr,       info)
      call prec%set('sub_restr',       p_choice%restr,      info)
      call prec%set('sub_prol',        p_choice%prol,       info)
      select case(trim(psb_toupper(p_choice%solve)))
      case('INVK')
        call prec%set('sub_solve',       p_choice%solve,   info)
      case('INVT')
        call prec%set('sub_solve',       p_choice%solve,   info)
      case('AINV')
        call prec%set('sub_solve',       p_choice%solve,   info)
        call prec%set('ainv_alg', p_choice%variant,   info)
      case default
        call prec%set('sub_solve',       p_choice%solve,   info)
        if (psb_toupper(p_choice%solve)=='MUMPS') &
             & call prec%set('mumps_loc_glob','local_solver',info)
      end select
      call prec%set('solver_sweeps',   p_choice%ssweeps,   info)
      call prec%set('sub_fillin',      p_choice%fill,       info)
      call prec%set('inv_fillin',      p_choice%invfill,    info)
      call prec%set('sub_iluthrs',     p_choice%thr,        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call prec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call prec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
      call prec%set('poly_degree',     p_choice%degree2,   info,pos='post')
      call prec%set('poly_variant',    p_choice%pvariant2, info,pos='post')
      if (p_choice%prhovalue > dzero ) then
         call prec%set('poly_rho_ba', p_choice%prhovalue2, info,pos='post')
      else
         call prec%set('poly_rho_estimate', p_choice%prhovariant2, info,pos='post')
      end if

      select case (psb_toupper(p_choice%smther2))
      case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI','L1-FBGS')
        ! do nothing
      case default
        call prec%set('sub_ovr',         p_choice%novr2,     info,pos='post')
        call prec%set('sub_restr',       p_choice%restr2,    info,pos='post')
        call prec%set('sub_prol',        p_choice%prol2,     info,pos='post')
        select case(trim(psb_toupper(p_choice%solve2)))
        case('INVK')
          call prec%set('sub_solve',       p_choice%solve2,   info)
        case('INVT')
          call prec%set('sub_solve',       p_choice%solve2,   info)
        case('AINV')
          call prec%set('sub_solve',       p_choice%solve2,   info)
          call prec%set('ainv_alg', p_choice%variant2,   info)
        case default
          call prec%set('sub_solve',       p_choice%solve2,   info, pos='post')
          if (psb_toupper(p_choice%solve2)=='MUMPS') &
               & call prec%set('mumps_loc_glob','local_solver',info)
        end select
        call prec%set('solver_sweeps',   p_choice%ssweeps2,   info,pos='post')
        call prec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call prec%set('inv_fillin',      p_choice%invfill2,  info,pos='post')
        call prec%set('sub_iluthrs',     p_choice%thr2,      info,pos='post')
      end select
    end if

    call prec%set('coarse_solve',    p_choice%csolve,    info)
    call prec%set('coarse_mat',      p_choice%cmat,      info)
     ! Set for the case of a KRM solver
    if (psb_toupper(p_choice%csolve) == 'KRM') then
      call prec%set('krm_method', p_choice%krm_method, info)
      call prec%set('krm_kprec', p_choice%krm_prec,   info)
      call prec%set('krm_sub_solve', p_choice%krm_subsolve, info)
      call prec%set('krm_global', p_choice%krm_global, info)
      call prec%set('krm_eps',   p_choice%krm_eps,   info)
      call prec%set('krm_irst',   p_choice%krm_irst,   info)
      call prec%set('krm_istopc',  p_choice%krm_istop,  info)
      call prec%set('krm_itmax',  p_choice%krm_itmax,  info)
      call prec%set('krm_itrace', p_choice%krm_itrace, info)
      call prec%set('krm_fillin', p_choice%krm_fillin, info)
    else
      if (psb_toupper(p_choice%csolve) == 'BJAC') &
        &  call prec%set('coarse_subsolve', p_choice%csbsolve,  info)
      call prec%set('coarse_fillin',   p_choice%cfill,     info)
      call prec%set('coarse_iluthrs',  p_choice%cthres,    info)
      call prec%set('coarse_sweeps',   p_choice%cjswp,     info)
    end if
  end select

  ! build the preconditioner
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  thier = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_hierarchy_bld')
    goto 9999
  end if
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info)
  tsmth = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ctxt, thier)
  call psb_amx(ctxt, tsmth)
  
  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')  thier+tsmth,thier,tsmth
    write(psb_out_unit,'(" ")')
  end if

  end if


  if (p_choice%dump) then
    if (mixed) then
      call sprec%dump(info,istart=p_choice%dlmin,iend=p_choice%dlmax,&
           & ac=p_choice%dump_ac,rp=p_choice%dump_rp,tprol=p_choice%dump_tprol,&
           & smoother=p_choice%dump_smoother, solver=p_choice%dump_solver, &
           & global_num=p_choice%dump_global_num)
    else
      call prec%dump(info,istart=p_choice%dlmin,iend=p_choice%dlmax,&
           & ac=p_choice%dump_ac,rp=p_choice%dump_rp,tprol=p_choice%dump_tprol,&
           & smoother=p_choice%dump_smoother, solver=p_choice%dump_solver, &
           & global_num=p_choice%dump_global_num)
    end if
  end if
  !
  ! iterative method parameters
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  select case(psb_toupper(trim(s_choice%kmethd)))
  case('RICHARDSON')
    call psb_richardson(a,prec,b,x,s_choice%eps,&
         & desc_a,info,itmax=s_choice%itmax,iter=iter,&
         & err=err,itrace=s_choice%itrace,&
         & istop=s_choice%istopc)
  case('BICGSTAB','BICGSTABL','BICG','CG','CGS','FCG','GCR','RGMRES')
    call psb_krylov(s_choice%kmethd,a,prec,b,x,s_choice%eps,&
         & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,itrace=s_choice%itrace,&
         & istop=s_choice%istopc,irst=s_choice%irst)
  case('CG-MIXED','CGMIXED')
    !
    ! Double precision CG with the single precision preconditioner
    ! built above on ASINGLE; see psb_cg_mixed_mod at the top of this file.
    !
    call psb_cg_mixed(a,sprec,b,x,s_choice%eps,&
         & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,&
         & itrace=s_choice%itrace,istop=s_choice%istopc)
  case default
    write(psb_err_unit,*) 'Unknown method :"',trim(s_choice%kmethd),'"'
    info=psb_err_invalid_input_
    call psb_errpush(info,name)
    goto 9999
  end select
  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1

  call psb_amx(ctxt,tslv)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1
  call psb_amx(ctxt,tslv)

  ! compute residual norms
  call psb_geall(r,desc_a,info)
  call r%zero()
  call psb_geasb(r,desc_a,info)
  call psb_geaxpby(done,b,dzero,r,desc_a,info)
  call psb_spmm(-done,a,x,done,r,desc_a,info)
  resmx  = psb_genrm2(r,desc_a,info)
  resmxp = psb_geamax(r,desc_a,info)

  vecsize = x%sizeof()
  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  if (mixed) then
    precsize = sprec%sizeof()
    nlv      = sprec%get_nlevs()
  else
    precsize = prec%sizeof()
    nlv      = prec%get_nlevs()
  end if
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,vecsize)
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)
  if (mixed) then
    call sprec%descr(info,iout=psb_out_unit)
  else
    call prec%descr(info,iout=psb_out_unit)
  end if
  if (iam == psb_root_) then
    write(psb_out_unit,'("Computed solution on ",i8," process(es)")')  np
    write(psb_out_unit,'("Number of threads                  : ",i12)') nth
    write(psb_out_unit,'("Total number of tasks              : ",i12)') nth*np
    write(psb_out_unit,'("Discretization domain size         : ",i12)') idim
    write(psb_out_unit,'("Linear system size                 : ",i12)') system_size
    write(psb_out_unit,'("PDE Coefficients                   : ",a)') trim(pdecoeff)
    write(psb_out_unit,'("Problem setup time                 : ",es12.5)') tpgen
    write(psb_out_unit,'("Krylov method                      : ",a)') trim(s_choice%kmethd)
    write(psb_out_unit,'("Preconditioner                     : ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Iterations to convergence          : ",i12)')    iter
    write(psb_out_unit,'("Relative error estimate on exit    : ",es12.5)') err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)')    nlv
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)') thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)') tsmth
    write(psb_out_unit,'("Total preconditioner setup time    : ",es12.5)') tsmth+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)') tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)') tslv/iter
    write(psb_out_unit,'("Total time                         : ",es12.5)') tslv+tsmth+thier
    write(psb_out_unit,'("Residual 2-norm                    : ",es12.5)') resmx
    write(psb_out_unit,'("Residual inf-norm                  : ",es12.5)') resmxp
    write(psb_out_unit,'("Total memory occupation for X      : ",i16)') vecsize
    write(psb_out_unit,'("Total memory occupation for A      : ",i16)') amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i16)') descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i16)') precsize
    write(psb_out_unit,'("Total memory occupation            : ",i16)') &
         & amatsize + descsize+precsize+2*vecsize    
    write(psb_out_unit,'("Storage format for A               : ",a  )') a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )') desc_a%get_fmt()

  end if
! call psb_print_timers(ctxt)
  !
  !  cleanup storage and exit
  !
  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  if (mixed) then
    call sprec%free(info)
    call psb_spfree(asingle,desc_a,info)
  else
    call prec%free(info)
  end if
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get iteration parameters from standard input
  !
  !
  ! get iteration parameters from standard input
  !
  subroutine get_parms(ctxt,afmt,idim,solve,prec,pdecoeff)

    implicit none

    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: idim
    character(len=*)    :: afmt
    type(solverdata)    :: solve
    type(precdata)      :: prec
    character(len=*)    :: pdecoeff
    integer(psb_ipk_)   :: iam, nm, np, inp_unit
    character(len=1024)   :: filename

    call psb_info(ctxt,iam,np)

    if (iam == psb_root_) then
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ctxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        write(psb_err_unit,*) 'Usage: amg_d_pde3d ctrl-file '
        call psb_abort(ctxt)
        stop        
      end if
      ! read input data
      !
      call read_data(afmt,inp_unit)            ! matrix storage format
      call read_data(idim,inp_unit)            ! Discretization grid size
      call read_data(pdecoeff,inp_unit)        ! PDE Coefficients
      ! Krylov solver data
      call read_data(solve%kmethd,inp_unit)    ! Krylov solver
      call read_data(solve%istopc,inp_unit)    ! stopping criterion
      call read_data(solve%itmax,inp_unit)     ! max num iterations
      call read_data(solve%itrace,inp_unit)    ! tracing
      call read_data(solve%irst,inp_unit)      ! restart
      call read_data(solve%eps,inp_unit)       ! tolerance
      ! preconditioner type
      call read_data(prec%descr,inp_unit)      ! verbose description of the prec
      call read_data(prec%ptype,inp_unit)      ! preconditioner type
      
      ! First smoother / 1-lev preconditioner
      call read_data(prec%smther,inp_unit)     ! smoother type
      call read_data(prec%jsweeps,inp_unit)    ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%degree,inp_unit)     ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%pvariant,inp_unit)   !
      call read_data(prec%prhovariant,inp_unit)! how to estimate rho(M^{-1}A)
      call read_data(prec%prhovalue,inp_unit)  ! if previous is set value, we set it from this one
      call read_data(prec%novr,inp_unit)       ! number of overlap layers
      call read_data(prec%restr,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve,inp_unit)      ! local subsolver
      call read_data(prec%ssweeps,inp_unit)    ! inner solver sweeps 
      call read_data(prec%variant,inp_unit)    ! AINV variant
      call read_data(prec%fill,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%invfill,inp_unit)    !Inverse fill-in for INVK
      call read_data(prec%thr,inp_unit)        ! threshold for ILUT
      ! Second smoother/ AMG post-smoother (if NONE ignored in main)
      call read_data(prec%smther2,inp_unit)     ! smoother type
      call read_data(prec%jsweeps2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%degree2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%pvariant2,inp_unit)   !
      call read_data(prec%prhovariant2,inp_unit)! how to estimate rho(M^{-1}A)
      call read_data(prec%prhovalue2,inp_unit)  ! if previous is set value, we set it from this one
      call read_data(prec%novr2,inp_unit)       ! number of overlap layers
      call read_data(prec%restr2,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol2,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve2,inp_unit)      ! local subsolver
      call read_data(prec%ssweeps2,inp_unit)    ! inner solver sweeps
      call read_data(prec%variant2,inp_unit)    ! AINV variant
      call read_data(prec%fill2,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%invfill2,inp_unit)    !Inverse fill-in for INVK
      call read_data(prec%thr2,inp_unit)        ! threshold for ILUT
      ! general AMG data
      call read_data(prec%mlcycle,inp_unit)     ! AMG cycle type
      call read_data(prec%outer_sweeps,inp_unit) ! number of 1lev/outer sweeps
      call read_data(prec%maxlevs,inp_unit)    ! max number of levels in AMG prec
      call read_data(prec%csizepp,inp_unit)       ! min size coarsest mat
      ! aggregation
      call read_data(prec%aggr_prol,inp_unit)    ! aggregation type
      call read_data(prec%par_aggr_alg,inp_unit)    ! parallel aggregation alg
      call read_data(prec%aggr_type,inp_unit)   ! type of aggregation
      call read_data(prec%aggr_size,inp_unit) ! Requested size of the aggregates for MATCHBOXP
      call read_data(prec%aggr_ord,inp_unit)    ! ordering for aggregation
      call read_data(prec%mncrratio,inp_unit)  ! minimum aggregation ratio
      call read_data(prec%aggr_filter,inp_unit) ! filtering
      call read_data(prec%athres,inp_unit)      ! smoothed aggr thresh
      call read_data(prec%thrvsz,inp_unit)      ! size of aggr thresh vector
      if (prec%thrvsz > 0) then
        call psb_realloc(prec%thrvsz,prec%athresv,info)
        call read_data(prec%athresv,inp_unit)   ! aggr thresh vector
      else
        read(inp_unit,*)                        ! dummy read to skip a record
      end if
      ! coasest-level solver
      call read_data(prec%csolve,inp_unit)      ! coarsest-lev solver
      call read_data(prec%csbsolve,inp_unit)    ! coarsest-lev subsolver
      call read_data(prec%cmat,inp_unit)        ! coarsest mat layout
      call read_data(prec%cfill,inp_unit)       ! fill-in for incompl LU
      call read_data(prec%cthres,inp_unit)      ! Threshold for ILUT
      call read_data(prec%cjswp,inp_unit)       ! sweeps for GS/JAC subsolver
      ! Krylov method for coarsest level
      call read_data(prec%krm_method,inp_unit)  ! Krylov method for coarsest level
      call read_data(prec%krm_prec,inp_unit)    ! Preconditioner for coarsest level
      call read_data(prec%krm_subsolve,inp_unit) ! Subsolver for coarsest level
      call read_data(prec%krm_global,inp_unit)  ! Is the solver global or local? TRUE or FALSE
      call read_data(prec%krm_eps,inp_unit)     ! Stopping tolerance
      call read_data(prec%krm_irst,inp_unit)    ! Restart for Krylov method
      call read_data(prec%krm_istop,inp_unit)   ! Stopping criterion
      call read_data(prec%krm_itmax,inp_unit)   ! Maximum number of iterations
      call read_data(prec%krm_itrace,inp_unit)  ! Trace of the lower Krylov iterations
      call read_data(prec%krm_fillin,inp_unit)  ! Fill-in for incomplete LU factorization
      ! dump 
      call read_data(prec%dump,inp_unit)       ! Dump on file?
      call read_data(prec%dlmin,inp_unit)      ! Minimum level to dump
      call read_data(prec%dlmax,inp_unit)      ! Maximum level to dump
      call read_data(prec%dump_ac,inp_unit)       
      call read_data(prec%dump_rp,inp_unit)       
      call read_data(prec%dump_tprol,inp_unit)    
      call read_data(prec%dump_smoother,inp_unit) 
      call read_data(prec%dump_solver,inp_unit)   
      call read_data(prec%dump_global_num,inp_unit)

      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    call psb_bcast(ctxt,afmt)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,pdecoeff)

    call psb_bcast(ctxt,solve%kmethd)
    call psb_bcast(ctxt,solve%istopc)
    call psb_bcast(ctxt,solve%itmax)
    call psb_bcast(ctxt,solve%itrace)
    call psb_bcast(ctxt,solve%irst)
    call psb_bcast(ctxt,solve%eps)

    call psb_bcast(ctxt,prec%descr)
    call psb_bcast(ctxt,prec%ptype)

    ! broadcast first (pre-)smoother / 1-lev prec data
    call psb_bcast(ctxt,prec%smther)
    call psb_bcast(ctxt,prec%jsweeps)
    call psb_bcast(ctxt,prec%degree)
    call psb_bcast(ctxt,prec%pvariant)
    call psb_bcast(ctxt,prec%prhovariant)
    call psb_bcast(ctxt,prec%prhovalue)  
    call psb_bcast(ctxt,prec%novr)
    call psb_bcast(ctxt,prec%restr)
    call psb_bcast(ctxt,prec%prol)
    call psb_bcast(ctxt,prec%solve)
    call psb_bcast(ctxt,prec%ssweeps)
    call psb_bcast(ctxt,prec%variant)
    call psb_bcast(ctxt,prec%fill)
    call psb_bcast(ctxt,prec%invfill)
    call psb_bcast(ctxt,prec%thr)
    ! broadcast second (post-)smoother
    call psb_bcast(ctxt,prec%smther2)
    call psb_bcast(ctxt,prec%jsweeps2)
    call psb_bcast(ctxt,prec%degree2)
    call psb_bcast(ctxt,prec%pvariant2)
    call psb_bcast(ctxt,prec%prhovariant2)
    call psb_bcast(ctxt,prec%prhovalue2)  
    call psb_bcast(ctxt,prec%novr2)
    call psb_bcast(ctxt,prec%restr2)
    call psb_bcast(ctxt,prec%prol2)
    call psb_bcast(ctxt,prec%solve2)
    call psb_bcast(ctxt,prec%ssweeps2)
    call psb_bcast(ctxt,prec%variant2)
    call psb_bcast(ctxt,prec%fill2)
    call psb_bcast(ctxt,prec%invfill2)
    call psb_bcast(ctxt,prec%thr2)

    ! broadcast AMG parameters
    call psb_bcast(ctxt,prec%mlcycle)
    call psb_bcast(ctxt,prec%outer_sweeps)
    call psb_bcast(ctxt,prec%maxlevs)
    call psb_bcast(ctxt,prec%csizepp)

    call psb_bcast(ctxt,prec%aggr_prol)
    call psb_bcast(ctxt,prec%par_aggr_alg)
    call psb_bcast(ctxt,prec%aggr_type)
    call psb_bcast(ctxt,prec%aggr_size)
    call psb_bcast(ctxt,prec%aggr_ord)
    call psb_bcast(ctxt,prec%aggr_filter)
    call psb_bcast(ctxt,prec%mncrratio)
    call psb_bcast(ctxt,prec%thrvsz)
    if (prec%thrvsz > 0) then
      if (iam /= psb_root_) call psb_realloc(prec%thrvsz,prec%athresv,info)
      call psb_bcast(ctxt,prec%athresv)
    end if
    call psb_bcast(ctxt,prec%athres)

    call psb_bcast(ctxt,prec%cmat)
    call psb_bcast(ctxt,prec%csolve)
    call psb_bcast(ctxt,prec%csbsolve)
    call psb_bcast(ctxt,prec%cfill)
    call psb_bcast(ctxt,prec%cthres)
    call psb_bcast(ctxt,prec%cjswp)
    ! Krylov method for coarsest level (broadcast)
    call psb_bcast(ctxt,prec%krm_method)
    call psb_bcast(ctxt,prec%krm_prec)
    call psb_bcast(ctxt,prec%krm_subsolve)
    call psb_bcast(ctxt,prec%krm_global)
    call psb_bcast(ctxt,prec%krm_eps)
    call psb_bcast(ctxt,prec%krm_irst)
    call psb_bcast(ctxt,prec%krm_istop)
    call psb_bcast(ctxt,prec%krm_itmax)
    call psb_bcast(ctxt,prec%krm_itrace)
    call psb_bcast(ctxt,prec%krm_fillin)
    ! dump
    call psb_bcast(ctxt,prec%dump)
    call psb_bcast(ctxt,prec%dlmin)
    call psb_bcast(ctxt,prec%dlmax)

    call psb_bcast(ctxt,prec%dump_ac)       
    call psb_bcast(ctxt,prec%dump_rp)       
    call psb_bcast(ctxt,prec%dump_tprol)    
    call psb_bcast(ctxt,prec%dump_smoother) 
    call psb_bcast(ctxt,prec%dump_solver)   
    call psb_bcast(ctxt,prec%dump_global_num)


    
  end subroutine get_parms

end program amg_d_pde3d

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
program amg_d_pde3d
  use psb_base_mod
  use amg_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use amg_d_pde3d_base_mod
  use amg_d_pde3d_exp_mod
  use amg_d_pde3d_gauss_mod
  use amg_d_genpde_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt, pdecoeff
  integer(psb_ipk_) :: idim
  integer(psb_epk_) :: system_size

  ! miscellaneous
  real(psb_dpk_) :: t1, t2, tprec, thier, tslv

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(amg_dprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: x,b,r
  ! parallel environment
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, resmx, resmxp

  ! Krylov solver data
  type solverdata
    character(len=40)  :: kmethd      ! Krylov solver
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
    character(len=16)  :: mlcycle      ! AMG cycle type
    integer(psb_ipk_)  :: maxlevs     ! maximum number of levels in AMG preconditioner

    ! AMG aggregation
    character(len=16)  :: aggr_prol    ! aggregation type: SMOOTHED, NONSMOOTHED
    character(len=16)  :: par_aggr_alg ! parallel aggregation algorithm: DEC, SYMDEC
    character(len=16)  :: aggr_type   ! Type of aggregation SOC1, SOC2, MATCHBOXP
    integer(psb_ipk_)  :: aggr_size   ! Requested size of the aggregates for MATCHBOXP
    character(len=16)  :: aggr_ord    ! ordering for aggregation: NATURAL, DEGREE
    character(len=16)  :: aggr_filter ! filtering: FILTER, NO_FILTER
    real(psb_dpk_)     :: mncrratio  ! minimum aggregation ratio
    integer(psb_ipk_)  :: matching_alg ! For NEW matching 1 2 3 variant
    real(psb_dpk_)     :: lambda       ! matching LAMBDA
    
    real(psb_dpk_), allocatable :: athresv(:) ! smoothed aggregation threshold vector
    integer(psb_ipk_)  :: thrvsz      ! size of threshold vector
    real(psb_dpk_)     :: athres      ! smoothed aggregation threshold
    integer(psb_ipk_)  :: csizepp     ! minimum size of coarsest matrix per process

    ! AMG smoother or pre-smoother; also 1-lev preconditioner
    character(len=16)  :: smther      ! (pre-)smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps     ! (pre-)smoother / 1-lev prec. sweeps
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    character(len=16)  :: restr       ! restriction over application of AS
    character(len=16)  :: prol        ! prolongation over application of AS
    character(len=16)  :: solve       ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    character(len=16)  :: variant     ! AINV variant: LLK, etc
    integer(psb_ipk_)  :: fill        ! fill-in for incomplete LU factorization
    integer(psb_ipk_)  :: invfill     ! Inverse fill-in for INVK
    real(psb_dpk_)     :: thr         ! threshold for ILUT factorization

    ! AMG post-smoother; ignored by 1-lev preconditioner
    character(len=16)  :: smther2     ! post-smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps2    ! post-smoother sweeps
    integer(psb_ipk_)  :: novr2       ! number of overlap layers
    character(len=16)  :: restr2      ! restriction  over application of AS
    character(len=16)  :: prol2       ! prolongation over application of AS
    character(len=16)  :: solve2      ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    character(len=16)  :: variant2    ! AINV variant: LLK, etc
    integer(psb_ipk_)  :: fill2       ! fill-in for incomplete LU factorization
    integer(psb_ipk_)  :: invfill2    ! Inverse fill-in for INVK
    real(psb_dpk_)      :: thr2        ! threshold for ILUT factorization

    ! coarsest-level solver
    character(len=16)  :: cmat        ! coarsest matrix layout: REPL, DIST
    character(len=16)  :: csolve      ! coarsest-lev solver: BJAC, SLUDIST (distr.
                                      ! mat.); UMF, MUMPS, SLU, ILU, ILUT, MILU
                                      ! (repl. mat.)
    character(len=16)  :: csbsolve    ! coarsest-lev local subsolver: ILU, ILUT,
                                      ! MILU, UMF, MUMPS, SLU
    integer(psb_ipk_)  :: cfill       ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: cthres      ! threshold for ILUT factorization
    integer(psb_ipk_)  :: cjswp       ! sweeps for GS or JAC coarsest-lev subsolver

  end type precdata
  type(precdata)       :: p_choice

  ! other variables
  integer(psb_ipk_)  :: info, i, k
  character(len=20)  :: name,ch_err

  info=psb_success_


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

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
  !  allocate and fill in the coefficient matrix, rhs and initial guess
  !

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  select case(psb_toupper(trim(pdecoeff)))
  case("CONST")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
        & a1,a2,a3,b1,b2,b3,c,g,info)
  case("EXP")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
        & a1_exp,a2_exp,a3_exp,b1_exp,b2_exp,b3_exp,c_exp,g_exp,info)
  case("GAUSS")
    call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
        & a1_gauss,a2_gauss,a3_gauss,b1_gauss,b2_gauss,b3_gauss,c_gauss,g_gauss,info)
  case default
    info=psb_err_from_subroutine_
    ch_err='amg_gen_pdecoeff'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select


  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='amg_gen_pde3d'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("PDE Coefficients             : ",a)')pdecoeff
  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')
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

  case ('BJAC')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case('AS')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_ovr',         p_choice%novr,    info)
    call prec%set('sub_restr',       p_choice%restr,   info)
    call prec%set('sub_prol',        p_choice%prol,    info)
    call prec%set('sub_solve',       p_choice%solve,   info)
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
    write(0,*) 'match variant ',p_choice%matching_alg
    if (p_choice%matching_alg>0)&
         & call prec%set('nwm_matching_alg',   p_choice%matching_alg, info)
    if (p_choice%lambda>0)&
         & call prec%set('nwm_lambda',   p_choice%lambda, info)

    call prec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call prec%set('aggr_filter',     p_choice%aggr_filter,info)


    call prec%set('smoother_type',   p_choice%smther,     info)
    call prec%set('smoother_sweeps', p_choice%jsweeps,    info)

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
      end select

      call prec%set('sub_fillin',      p_choice%fill,       info)
      call prec%set('inv_fillin',      p_choice%invfill,    info)
      call prec%set('sub_iluthrs',     p_choice%thr,        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call prec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call prec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
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
        end select

        call prec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call prec%set('inv_fillin',      p_choice%invfill2,  info,pos='post')
        call prec%set('sub_iluthrs',     p_choice%thr2,      info,pos='post')
      end select
    end if

    call prec%set('coarse_solve',    p_choice%csolve,    info)
    if (psb_toupper(p_choice%csolve) == 'BJAC') &
         &  call prec%set('coarse_subsolve', p_choice%csbsolve,  info)
    call prec%set('coarse_mat',      p_choice%cmat,      info)
    call prec%set('coarse_fillin',   p_choice%cfill,     info)
    call prec%set('coarse_iluthrs',  p_choice%cthres,    info)
    call prec%set('coarse_sweeps',   p_choice%cjswp,     info)

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
  tprec = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ctxt, thier)
  call psb_amx(ctxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')thier+tprec
    write(psb_out_unit,'(" ")')
  end if

  !
  ! iterative method parameters
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_krylov(s_choice%kmethd,a,prec,b,x,s_choice%eps,&
       & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,itrace=s_choice%itrace,&
       & istop=s_choice%istopc,irst=s_choice%irst)
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

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)
  call prec%descr(info,iout=psb_out_unit)
  if (iam == psb_root_) then
    write(psb_out_unit,'("Computed solution on ",i8," processors")')  np
    write(psb_out_unit,'("Linear system size                 : ",i12)') system_size
    write(psb_out_unit,'("PDE Coefficients                   : ",a)') trim(pdecoeff)
    write(psb_out_unit,'("Krylov method                      : ",a)') trim(s_choice%kmethd)
    write(psb_out_unit,'("Preconditioner                     : ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Iterations to convergence          : ",i12)')    iter
    write(psb_out_unit,'("Relative error estimate on exit    : ",es12.5)') err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)')    prec%get_nlevs()
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)') thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)') tprec
    write(psb_out_unit,'("Total time for preconditioner      : ",es12.5)') tprec+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)') tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)') tslv/iter
    write(psb_out_unit,'("Total time                         : ",es12.5)') tslv+tprec+thier
    write(psb_out_unit,'("Residual 2-norm                    : ",es12.5)') resmx
    write(psb_out_unit,'("Residual inf-norm                  : ",es12.5)') resmxp
    write(psb_out_unit,'("Total memory occupation for A      : ",i12)') amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i12)') descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i12)') precsize
    write(psb_out_unit,'("Storage format for A               : ",a  )') a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )') desc_a%get_fmt()

  end if

  !
  !  cleanup storage and exit
  !
  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call prec%free(info)
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
        inp_unit=psb_inp_unit
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
      call read_data(prec%novr,inp_unit)       ! number of overlap layers
      call read_data(prec%restr,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve,inp_unit)      ! local subsolver
      call read_data(prec%variant,inp_unit)    ! AINV variant
      call read_data(prec%fill,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%invfill,inp_unit)    !Inverse fill-in for INVK
      call read_data(prec%thr,inp_unit)        ! threshold for ILUT
      ! Second smoother/ AMG post-smoother (if NONE ignored in main)
      call read_data(prec%smther2,inp_unit)     ! smoother type
      call read_data(prec%jsweeps2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%novr2,inp_unit)       ! number of overlap layers
      call read_data(prec%restr2,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol2,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve2,inp_unit)      ! local subsolver
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
      call read_data(prec%aggr_filter,inp_unit) ! filtering
      call read_data(prec%mncrratio,inp_unit)  ! minimum aggregation ratio
      call read_data(prec%matching_alg,inp_unit)  ! matching variant
      call read_data(prec%lambda,inp_unit)      ! lambda 
      call read_data(prec%thrvsz,inp_unit)      ! size of aggr thresh vector
      if (prec%thrvsz > 0) then
        call psb_realloc(prec%thrvsz,prec%athresv,info)
        call read_data(prec%athresv,inp_unit)   ! aggr thresh vector
      else
        read(inp_unit,*)                        ! dummy read to skip a record
      end if
      call read_data(prec%athres,inp_unit)      ! smoothed aggr thresh
      ! coasest-level solver
      call read_data(prec%csolve,inp_unit)      ! coarsest-lev solver
      call read_data(prec%csbsolve,inp_unit)    ! coarsest-lev subsolver
      call read_data(prec%cmat,inp_unit)        ! coarsest mat layout
      call read_data(prec%cfill,inp_unit)       ! fill-in for incompl LU
      call read_data(prec%cthres,inp_unit)      ! Threshold for ILUT
      call read_data(prec%cjswp,inp_unit)       ! sweeps for GS/JAC subsolver
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
    call psb_bcast(ctxt,prec%novr)
    call psb_bcast(ctxt,prec%restr)
    call psb_bcast(ctxt,prec%prol)
    call psb_bcast(ctxt,prec%solve)
    call psb_bcast(ctxt,prec%variant)
    call psb_bcast(ctxt,prec%fill)
    call psb_bcast(ctxt,prec%invfill)
    call psb_bcast(ctxt,prec%thr)
    ! broadcast second (post-)smoother
    call psb_bcast(ctxt,prec%smther2)
    call psb_bcast(ctxt,prec%jsweeps2)
    call psb_bcast(ctxt,prec%novr2)
    call psb_bcast(ctxt,prec%restr2)
    call psb_bcast(ctxt,prec%prol2)
    call psb_bcast(ctxt,prec%solve2)
    call psb_bcast(ctxt,prec%variant2)
    call psb_bcast(ctxt,prec%fill2)
    call psb_bcast(ctxt,prec%invfill2)
    call psb_bcast(ctxt,prec%thr2)

    ! broadcast AMG parameters
    call psb_bcast(ctxt,prec%mlcycle)
    call psb_bcast(ctxt,prec%outer_sweeps)
    call psb_bcast(ctxt,prec%maxlevs)

    call psb_bcast(ctxt,prec%aggr_prol)
    call psb_bcast(ctxt,prec%par_aggr_alg)
    call psb_bcast(ctxt,prec%aggr_type)
    call psb_bcast(ctxt,prec%aggr_size)
    call psb_bcast(ctxt,prec%aggr_ord)
    call psb_bcast(ctxt,prec%aggr_filter)
    call psb_bcast(ctxt,prec%mncrratio)
    call psb_bcast(ctxt,prec%matching_alg)
    call psb_bcast(ctxt,prec%lambda)


    call psb_bcast(ctxt,prec%thrvsz)
    if (prec%thrvsz > 0) then
      if (iam /= psb_root_) call psb_realloc(prec%thrvsz,prec%athresv,info)
      call psb_bcast(ctxt,prec%athresv)
    end if
    call psb_bcast(ctxt,prec%athres)

    call psb_bcast(ctxt,prec%csizepp)
    call psb_bcast(ctxt,prec%cmat)
    call psb_bcast(ctxt,prec%csolve)
    call psb_bcast(ctxt,prec%csbsolve)
    call psb_bcast(ctxt,prec%cfill)
    call psb_bcast(ctxt,prec%cthres)
    call psb_bcast(ctxt,prec%cjswp)


  end subroutine get_parms

end program amg_d_pde3d

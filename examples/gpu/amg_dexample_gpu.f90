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
! File: amg_dexample_gpu.f90
!
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. The solver is CG, coupled with one of the
! following multi-level preconditioner, as explained in Section 4.2 of
! the AMG4PSBLAS User's and Reference Guide:
!
! - choice = 1, a V-cycle with decoupled smoothed aggregation, 4 Jacobi
! sweeps as pre/post-smoother and 8 Jacobi sweeps as coarsest-level
! solver with replicated coarsest matrix
!
! - choice = 2, a W-cycle based on the coupled aggregation relying on matching,
! with maximum size of aggregates equal to 8 and smoothed prolongators,
! 2 sweeps of Block-Jacobi ipre/post-smoother using approximate inverse INVK and
! 4 sweeps of Block-Jacobi with INVK as coarsest-level solver on distributed
! coarsest matrix 
!
! The matrix and the rhs are read from files (if an rhs is not available, the
! unit rhs is set).
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
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
program amg_dexample_gpu
  use psb_base_mod
  use amg_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use psb_gpu_mod
  use data_input
  use amg_d_pde_mod
  implicit none

  ! input parameters

  ! sparse matrices
  type(psb_dspmat_type) :: A

  ! sparse matrices descriptor
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(amg_dprec_type)  :: P

  ! right-hand side, solution and residual vectors
  type(psb_d_vect_type) :: x, b, r
  ! GPU variables
  type(psb_d_hlg_sparse_mat) :: agmold
  type(psb_d_vect_gpu)       :: vgmold
  type(psb_i_vect_gpu)       :: igmold

  ! solver and preconditioner parameters
  real(psb_dpk_)   :: tol, err
  integer          :: itmax, iter, istop
  integer          :: nlev

  ! parallel environment parameters
  type(psb_ctxt_type) :: ctxt
  integer             :: iam, np

  ! other variables
  integer            :: choice       
  integer            :: i,info,j
  integer(psb_epk_) :: amatsize, precsize, descsize
  integer(psb_epk_) :: system_size
  integer            :: idim, ierr, ircode
  real(psb_dpk_)     :: resmx, resmxp
  real(psb_dpk_)     :: t1, t2, tprec
  character(len=5)   :: afmt='CSR'
  character(len=20)  :: name, kmethod

  ! initialize the parallel environment

  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
  !
  ! BEWARE: if you have NGPUS  per node, the default is to
  ! attach to mod(IAM,NGPUS)
  !
  call psb_gpu_init(ictxt)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif

  name='amg_dexample_gpu'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to AMG4PSBLAS version: ',amg_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
  write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
  write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_gpu_DeviceName())  

  ! get parameters

  call get_parms(ctxt,choice,idim,itmax,tol)

  !  allocate and fill in the coefficient matrix, rhs and initial guess

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info)  
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(*,'(" ")')

  select case(choice)

 case(1)

    ! initialize a V-cycle preconditioner, relying on decoupled smoothed aggregation 
    ! with 4 Jacobi sweeps as pre/post-smoother
    ! and 8 Jacobi sweeps as coarsest-level solver on replicated coarsest matrix

    call P%init(ctxt,'ML',info)
    call P%set('SMOOTHER_TYPE','JACOBI',info)
    call P%set('SMOOTHER_SWEEPS',4,info)
    call P%set('COARSE_SOLVE','JACOBI',info)
    call P%set('COARSE_SWEEPS',8,info)
    kmethod = 'CG'

  case(2)

   ! initialize a W-cycle preconditioner based on the coupled aggregation relying on matching,
   ! with maximum size of aggregates equal to 8 and smoothed prolongators,
   ! 2 sweeps of Block-Jacobi pre/post-smoother using approximate inverse INVK and
   ! 4 sweeps of Block-Jacobi with INVK on the coarsest level distributed matrix

    call P%init(ctxt,'ML',info)
    call P%set('PAR_AGGR_ALG','COUPLED',info)
    call P%set('AGGR_TYPE','MATCHBOXP',info)
    call P%set('AGGR_SIZE',8,info)
    call P%set('ML_CYCLE','WCYCLE',info)
    call P%set('SMOOTHER_TYPE','BJAC',info)
    call P%set('SMOOTHER_SWEEPS',2,info)
    call P%set('SUB_SOLVE','INVK',info)  
    call P%set('COARSE_SOLVE','BJAC',info)
    call P%set('COARSE_SUBSOLVE','INVK',info)
    call P%set('COARSE_SWEEPS',4,info)
    call P%set('COARSE_MAT','DIST',info)
    kmethod = 'CG'

  end select

  call psb_barrier(ctxt)
  t1 = psb_wtime()

  ! build the preconditioner
  call P%hierarchy_build(A,desc_A,info)
  call P%smoothers_build(A,desc_A,info, amold=agmold, vmold=vgmold, imold=igmold)

  tprec = psb_wtime()-t1
  call psb_amx(ctxt, tprec)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_precbld')
    goto 9999
  end if

  ! set the solver parameters and the initial guess

  call psb_geall(x,desc_A,info)
  call x%zero()
  call psb_geasb(x,desc_A,info)

  ! Convert A, DESC_A,X,B to a GPU-enabled format
  call desc_a%cnv(mold=igmold)
  call a%cscnv(info,mold=agmold)
  call psb_geasb(x,desc_a,info,mold=vgmold)
  call psb_geasb(b,desc_a,info,mold=vgmold)


  ! solve Ax=b with preconditioned Krylov method

  call psb_barrier(ctxt)
  call prec%allocate_wrk(info)
  t1 = psb_wtime()

  call psb_krylov(kmethod,A,P,b,x,tol,desc_A,info,itmax,iter,err,itrace=1,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)
  call prec%deallocate_wrk(info)

  call psb_geall(r,desc_A,info)
  call r%zero()
  call psb_geasb(r,desc_A,info)
  call psb_geaxpby(done,b,dzero,r,desc_A,info)
  call psb_spmm(-done,A,x,done,r,desc_A,info)
  resmx  = psb_genrm2(r,desc_A,info)
  resmxp = psb_geamax(r,desc_A,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = p%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)

  call P%descr(info)

  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Matrix from PDE example")')
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Linear system size        : ",i12)') system_size
    write(*,'("Krylov method             : ",a)') kmethod
    write(*,'("Iterations to convergence : ",i6)')iter
    write(*,'("Error estimate on exit    : ",es12.5)')err
    write(*,'("Time to build prec.       : ",es12.5)')tprec
    write(*,'("Time to solve system      : ",es12.5)')t2
    write(*,'("Time per iteration        : ",es12.5)')t2/(iter)
    write(*,'("Total time                : ",es12.5)')t2+tprec
    write(*,'("Residual 2-norm           : ",es12.5)')resmx
    write(*,'("Residual inf-norm         : ",es12.5)')resmxp
    write(*,'("Total memory occupation for A      : ",i12)')amatsize
    write(*,'("Total memory occupation for DESC_A : ",i12)')descsize
    write(*,'("Total memory occupation for PREC   : ",i12)')precsize
  end if

  call psb_gefree(b, desc_A,info)
  call psb_gefree(x, desc_A,info)
  call psb_spfree(A, desc_A,info)
  call P%free(info)
  call psb_cdfree(desc_A,info)
  call psb_gpu_exit()
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ctxt,choice,idim,itmax,tol)

    implicit none

    type(psb_ctxt_type) :: ctxt
    integer             :: choice, idim, itmax
    real(psb_dpk_)      :: tol
    integer             :: iam, np, inp_unit
    character(len=1024) :: filename

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
      ! read input parameters
      call read_data(choice,inp_unit)
      call read_data(idim,inp_unit)
      call read_data(itmax,inp_unit)
      call read_data(tol,inp_unit)
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if   
    end if

    call psb_bcast(ctxt,choice)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,itmax)
    call psb_bcast(ctxt,tol)

  end subroutine get_parms

end program amg_dexample_gpu

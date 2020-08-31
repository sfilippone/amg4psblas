
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
! File: amg_dexample_1lev.f90
!
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. The solver is BiCGStab preconditioned by
! RAS with overlap 2 and ILU(0) on the local blocks, as explained in Section 5.1 
! of the MLD2P4 User's and Reference Guide.
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
program amg_dexample_1lev
  use psb_base_mod
  use amg_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use amg_d_pde_mod
  implicit none


  ! sparse matrices
  type(psb_dspmat_type) :: A

  ! descriptor of sparse matrices
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(amg_dprec_type)  :: P

  ! right-hand side, solution and residual vectors
  type(psb_d_vect_type) :: x, b, r

  ! solver parameters
  real(psb_dpk_)   :: tol, err
  integer :: itmax, iter, itrace, istop

  ! parallel environment parameters
  integer            :: ictxt, iam, np

  ! other variables
  integer            :: i,info,j
  integer(psb_epk_) :: amatsize, precsize, descsize
  integer(psb_epk_) :: system_size
  integer            :: idim, nlev, ierr, ircode
  real(psb_dpk_) :: resmx, resmxp
  real(psb_dpk_) :: t1, t2, tprec
  character(len=5)   :: afmt='CSR'
  character(len=20)  :: name, kmethod

  ! initialize the parallel environment
  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif

  name='amg_dexample_ml'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to MLD2P4 version: ',amg_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  ! get parameters

  call get_parms(ictxt,idim,itmax,tol)

  ! allocate and fill in the coefficient matrix, rhs and initial guess

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call amg_gen_pde3d(ictxt,idim,a,b,x,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(*,'(" ")')

  ! set RAS

  call P%init(ictxt,'AS',info)

  ! set number of overlaps

  call P%set('SUB_OVR',2,info)

  ! build the preconditioner

  t1 = psb_wtime()

  call P%build(A,desc_A,info)

  tprec = psb_wtime()-t1
  call psb_amx(ictxt, tprec)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_precbld')
    goto 9999
  end if

  ! set the initial guess

  call psb_geall(x,desc_A,info)
  call x%zero()
  call psb_geasb(x,desc_A,info)

  ! solve Ax=b with preconditioned Krylov method: BiCGSTAB
  kmethod = 'BiCGSTAB'

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  call psb_krylov(kmethod,A,P,b,x,tol,desc_A,info,itmax,iter,err,itrace=1,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)

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
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)

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
  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ictxt,idim,itmax,tol)

    implicit none

    integer             :: idim, ictxt, itmax
    real(psb_dpk_)      :: tol
    integer             :: iam, np, inp_unit
    character(len=1024)   :: filename

    call psb_info(ictxt,iam,np)
    
    if (iam == psb_root_) then
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ictxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      
      ! read input parameters
      call read_data(idim,inp_unit)
      call read_data(itmax,inp_unit)
      call read_data(tol,inp_unit)
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,tol)

  end subroutine get_parms
end program amg_dexample_1lev

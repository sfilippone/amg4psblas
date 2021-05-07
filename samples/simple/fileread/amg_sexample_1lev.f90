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
! File: amg_sexample_1lev.f90
!
! This sample program solves a linear system by using BiCGStab preconditioned by
! RAS with overlap 2 and ILU(0) on the local blocks, as explained in Section 4.1 
! of the AMG4PSBLAS User's and Reference Guide.
!
! The matrix and the rhs are read from files (if an rhs is not available, the
! unit rhs is set). 
!
program amg_sexample_1lev
  use psb_base_mod
  use amg_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  implicit none

  ! input parameters
  character(len=40) :: mtrx_file, rhs_file
  character(len=2)  :: filefmt

  ! sparse matrices
  type(psb_sspmat_type) :: A, aux_A

  ! descriptor of sparse matrices
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(amg_sprec_type)  :: P

  ! right-hand side, solution and residual vectors
  type(psb_s_vect_type)  :: b, x, r
  real(psb_spk_), allocatable , save  :: x_glob(:), r_glob(:)
  real(psb_spk_), allocatable, target :: aux_b(:,:)
  real(psb_spk_), pointer  :: b_glob(:)

  ! solver and preconditioner parameters
  real(psb_spk_)   :: tol, err
  integer          :: itmax, iter, istop
  integer          :: nlev

  ! parallel environment parameters
  type(psb_ctxt_type) :: ctxt
  integer             :: iam, np

  ! other variables
  integer            :: i,info,j,m_problem
  integer(psb_epk_) :: amatsize, precsize, descsize
  integer :: ierr, ircode
  real(psb_spk_) :: resmx, resmxp
  real(psb_dpk_) :: t1, t2, tprec
  character(len=20)  :: name, kmethod
  integer, parameter :: iunit=12

  ! initialize the parallel environment

  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif

  name='amg_sexample_ml'
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

  !  get parameters

  call get_parms(ctxt,mtrx_file,rhs_file,filefmt,itmax,tol)

  call psb_barrier(ctxt)
  t1 = psb_wtime()  

  ! read and assemble the matrix A and the right-hand side b
  ! using PSBLAS routines for sparse matrix / vector management

  if (iam == psb_root_) then
    select case(psb_toupper(filefmt)) 
    case('MM') 
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS. 
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if (info == psb_success_) then 
        if (rhs_file /= 'NONE') then
          call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
        end if
      end if

    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)

    case default
      info = -1 
      write(0,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= psb_success_) then
      write(0,*) 'Error while reading input matrix '
      call psb_abort(ctxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ctxt,m_problem)

    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,1) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(0,'("Ok, got an rhs ")')
      b_glob =>aux_b(:,1)
    else
      write(*,'("Generating an rhs...")')
      write(*,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      endif

      b_glob => aux_b(:,1)
      do i=1, m_problem
        b_glob(i) = 1.d0
      enddo
    endif
  else
    call psb_bcast(ctxt,m_problem)
  end if

  call psb_barrier(ctxt)
  if (iam == psb_root_) write(*,'("Partition type: block")')
  call psb_matdist(aux_A, A, ctxt, desc_A,info,parts=part_block)
  call psb_scatter(b_glob,b,desc_a,info,root=psb_root_)

  t2 = psb_wtime() - t1

  call psb_amx(ctxt, t2)

  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to read and partition matrix : ",es12.5)')t2
    write(*,'(" ")')
  end if


  ! START SETTING PARAMETER

  ! set RAS

  call P%init(ctxt,'AS',info)

  ! set number of overlaps

  call P%set('SUB_OVR',2,info)
    
  ! build the preconditioner

  t1 = psb_wtime()

  call P%build(A,desc_A,info)

  tprec = psb_wtime()-t1
  call psb_amx(ctxt, tprec)

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

  call psb_barrier(ctxt)
  t1 = psb_wtime()

  call psb_krylov(kmethod,A,P,b,x,tol,desc_A,info,itmax,iter,err,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)

  call psb_geasb(r,desc_A,info,scratch=.true.)
  call psb_geaxpby(sone,b,szero,r,desc_A,info)
  call psb_spmm(-sone,A,x,sone,r,desc_A,info)
  resmx  = psb_genrm2(r,desc_A,info)
  resmxp = psb_geamax(r,desc_A,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = p%sizeof()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)

  call P%descr(info)
  if (iam == psb_root_) then 
    write(*,'(" ")')
    write(*,'("Matrix: ",A)')mtrx_file
    write(*,'("Computed solution on ",i8," processors")')np
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

  call psb_gather(x_glob,x,desc_a,info,root=psb_root_)
  if (info == psb_success_) &
       & call psb_gather(r_glob,r,desc_a,info,root=psb_root_)
  if (info /= psb_success_) goto 9999
  if (iam == psb_root_) then
    write(0,'(" ")')
    write(0,'("Saving x on file")')
    write(20,*) 'matrix: ',mtrx_file
    write(20,*) 'computed solution on ',np,' processor(s).'
    write(20,*) 'iterations to convergence: ',iter
    write(20,*) 'error estimate (infinity norm) on exit:', &
         & ' ||r||/(||a||||x||+||b||) = ',err
    write(20,*) 'max residual = ',resmx, resmxp
    write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
    do i=1,m_problem
      write(20,998) i,x_glob(i),r_glob(i),b_glob(i)
    enddo
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  ! deallocate the data structures

  call psb_gefree(b, desc_A,info)
  call psb_gefree(x, desc_A,info)
  call psb_spfree(A, desc_A,info)
  call P%free(info)
  call psb_cdfree(desc_A,info)
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ctxt,mtrx,rhs,filefmt,itmax,tol)

    implicit none

    type(psb_ctxt_type) :: ctxt
    integer             :: itmax
    real(psb_spk_)      :: tol
    character(len=*)    :: mtrx, rhs,filefmt
    integer             :: iam, np, inp_unit
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
      ! read input parameters
      call read_data(mtrx,inp_unit)
      call read_data(rhs,inp_unit)
      call read_data(filefmt,inp_unit)
      call read_data(itmax,inp_unit)
      call read_data(tol,inp_unit)
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if     
    end if

    call psb_bcast(ctxt,mtrx)
    call psb_bcast(ctxt,rhs)
    call psb_bcast(ctxt,filefmt)
    call psb_bcast(ctxt,itmax)
    call psb_bcast(ctxt,tol)

  end subroutine get_parms
end program amg_sexample_1lev

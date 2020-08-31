!   
!   
!                             MLD2P4  version 2.2
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
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
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
! File: mld_d_pde3d.f90
!
! Program: mld_d_pde3d
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
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
module mld_d_pde3d_mod
contains
  !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=dzero/sqrt((3*done))
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=dzero/sqrt((3*done))
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      
    b3=dzero/sqrt((3*done))
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z      
    c=dzero
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  a1   
    real(psb_dpk_), intent(in) :: x,y,z
    a1=done!/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=done!/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=done!/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_,done,dzero
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then 
      g = exp(y**2-z**2)
    end if
  end function g
end module mld_d_pde3d_mod

program mld_d_pde3d
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use mld_d_pde3d_mod
  use mld_d_tlu_solver
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim

  ! miscellaneous 
  real(psb_dpk_) :: t1, t2, tprec, thier, tslv

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(mld_dprec_type)  :: prec
  type(mld_d_tlu_solver_type) :: tlusv
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: x,b
  ! parallel environment
  integer(psb_ipk_) :: ictxt, iam, np

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps

  ! other variables
  integer(psb_ipk_)  :: info, i
  character(len=20)  :: name,ch_err

  info=psb_success_


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='mld_d_pde3d'
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to MLD2P4 version: ',mld_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  !
  !  get parameters
  !
  call get_parms(ictxt,kmethd,afmt,idim,istopc,itmax,itrace,irst,eps)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_gen_pde3d(ictxt,idim,a,b,x,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')
  !
  !  prepare the preconditioner: an ML with defaults, but with TLU solver at
  !  intermediate levels. All other parameters are at default values. 
  !  
  call prec%init('ML',       info)
  
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='prec%hierarchy_bld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  thier = psb_wtime()-t1
  nlv = prec%get_nlevs()
  call prec%set(tlusv,   info,ilev=1,ilmax=max(1,nlv-1))
  
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='prec%smoothers_build'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  tprec = psb_wtime()-t1

  call psb_amx(ictxt,thier)
  call psb_amx(ictxt,tprec)

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Preconditioner time : ",es12.5)') tprec+thier
  if (iam == psb_root_) call prec%descr(info)
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')

  !
  ! iterative method parameters 
  !
  if(iam == psb_root_) &
       & write(psb_out_unit,'("Calling iterative method ",a)')kmethd
  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  call psb_krylov(kmethd,a,prec,b,x,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)     

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ictxt)
  tslv = psb_wtime() - t1
  call psb_amx(ictxt,tslv)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Numer of levels of aggr. hierarchy: ",i12)') prec%get_nlevs()
    write(psb_out_unit,'("Time to build aggr. hierarchy     : ",es12.5)')  thier
    write(psb_out_unit,'("Time to build smoothers           : ",es12.5)')  tprec
    write(psb_out_unit,'("Total preconditioner time         : ",es12.5)')  tprec+thier
    write(psb_out_unit,'("Time to solve system              : ",es12.5)')  tslv
    write(psb_out_unit,'("Time per iteration                : ",es12.5)')  tslv/iter
    write(psb_out_unit,'("Number of iterations              : ",i0)')      iter
    write(psb_out_unit,'("Convergence indicator on exit     : ",es12.5)')  err
    write(psb_out_unit,'("Info  on exit                     : ",i0)')      info
    write(psb_out_unit,'("Total memory occupation for      A: ",i12)') amatsize
    write(psb_out_unit,'("Storage format for               A: ",a)')   trim(a%get_fmt())
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)') descsize
    write(psb_out_unit,'("Storage format for          DESC_A: ",a)')   trim(desc_a%get_fmt())
    write(psb_out_unit,'("Total memory occupation for   PREC: ",i12)') precsize
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
  call psb_exit(ictxt)
  stop

9999 continue
  call psb_error(ictxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,kmethd,afmt,idim,istopc,itmax,itrace,irst,eps)
    
    integer(psb_ipk_) :: ictxt
    character(len=*)  :: kmethd, afmt
    integer(psb_ipk_) :: idim, istopc,itmax,itrace,irst
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: eps
    character(len=20) :: buffer

    call psb_info(ictxt, iam, np)

    if (iam == psb_root_) then
      call read_data(kmethd,psb_inp_unit)
      call read_data(afmt,psb_inp_unit)
      call read_data(idim,psb_inp_unit)
      call read_data(istopc,psb_inp_unit)
      call read_data(itmax,psb_inp_unit)
      call read_data(itrace,psb_inp_unit)
      call read_data(irst,psb_inp_unit)
      call read_data(eps,psb_inp_unit)
    end if

    ! broadcast parameters to all processors
    call psb_bcast(ictxt,kmethd)
    call psb_bcast(ictxt,afmt)
    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,istopc)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,itrace)
    call psb_bcast(ictxt,irst)
    call psb_bcast(ictxt,eps)

    if (iam == psb_root_) then 
      write(psb_out_unit,'("Solving matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)') np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'("Preconditioner       : ",a)') 'ML-TLU'
      write(psb_out_unit,'("Iterative method     : ",a)') kmethd
      write(psb_out_unit,'(" ")')
    endif

    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer(psb_ipk_) :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  mld_d_pde3d methd prec dim &
         &[istop itmax itrace]'  
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl' 
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '  
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage

end program mld_d_pde3d

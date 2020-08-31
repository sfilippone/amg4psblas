!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.0)
!    
!    (C) Copyright 2008,2009,2010,2012,2013
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
!  Current version of this file contributed by:
!        Ambra Abdullahi Hassan 
!
!  
subroutine d_mumps_solver_apply(alpha,sv,x,beta,y,desc_data,&
     & trans,work,info,init,initu)
  use psb_base_mod
  use amg_d_mumps_solver
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(amg_d_mumps_solver_type), intent(inout) :: sv
  real(psb_dpk_),intent(inout)         :: x(:)
  real(psb_dpk_),intent(inout)         :: y(:)
  real(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)          :: trans
  real(psb_dpk_),target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(out)       :: info
  character, intent(in), optional       :: init
  real(psb_dpk_),intent(inout), optional :: initu(:)

  integer(psb_ipk_)  :: n_row, n_col
  integer(psb_lpk_)  :: nglob
  integer(psb_epk_)  :: eng
  real(psb_dpk_), allocatable     :: ww(:)
  real(psb_dpk_), allocatable, target :: gx(:)
  integer(psb_ipk_)  :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_mumps_solver_apply'

  call psb_erractionsave(err_act)

#if defined(HAVE_MUMPS_) 
  info = psb_success_
  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select
  !
  ! For non-iterative solvers, init and initu are ignored.
  !

  nglob = desc_data%get_global_rows()
  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  ! Running in local mode? 
  if (sv%ipar(1) == amg_local_solver_ ) then
    gx = x
  else if (sv%ipar(1) == amg_global_solver_ ) then

    if (n_col <= size(work)) then 
      ww = work(1:n_col)
    else
      allocate(ww(n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/n_col/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if
    end if
    allocate(gx(nglob),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_; eng = nglob
      call psb_errpush(info,name,e_err=(/eng/),&
           & a_err='real(psb_dpk_)')
      goto 9999      
    end if
    call psb_gather(gx, x, desc_data, info, root=izero)
  else
      info=psb_err_internal_error_
      call psb_errpush(info,name,& 
           & a_err='Invalid local/global solver in MUMPS')
      goto 9999      
  end if
  
  select case(trans_)
  case('N')
    sv%id%icntl(9) = 1
  case('T')
    sv%id%icntl(9) = 2
  case default
    call psb_errpush(psb_err_internal_error_,&
         & name,a_err='Invalid TRANS in subsolve')
    goto 9999
  end select

  sv%id%rhs  => gx
  sv%id%nrhs =  1
  sv%id%icntl(1)=-1
  sv%id%icntl(2)=-1
  sv%id%icntl(3)=-1
  sv%id%icntl(4)=-1
  sv%id%job = 3
  call dmumps(sv%id)
  if (sv%ipar(1) == amg_local_solver_ ) then
    call psb_geaxpby(alpha,gx,beta,y,desc_data,info)
  else
    call psb_scatter(gx, ww, desc_data, info, root=izero)   
    if (info == psb_success_) then
      call psb_geaxpby(alpha,ww,beta,y,desc_data,info)
    end if
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,&
         & name,a_err='Error in subsolve')
    goto 9999
  endif

  if (allocated(ww)) deallocate(ww)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

#else
  write(psb_err_unit,*) "MUMPS Not Configured, fix make.inc and recompile "
#endif
end subroutine d_mumps_solver_apply


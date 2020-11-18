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
subroutine amg_s_diag_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
     & trans,work,wv,info,init,initu)
  
  use psb_base_mod
  use amg_s_diag_solver, amg_protect_name => amg_s_diag_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)                :: desc_data
  class(amg_s_diag_solver_type), intent(inout) :: sv
  type(psb_s_vect_type), intent(inout)         :: x
  type(psb_s_vect_type), intent(inout)         :: y
  real(psb_spk_),intent(in)                     :: alpha,beta
  character(len=1),intent(in)                    :: trans
  real(psb_spk_),target, intent(inout)          :: work(:)
  type(psb_s_vect_type),intent(inout)          :: wv(:)
  integer(psb_ipk_), intent(out)                 :: info
  character, intent(in), optional                :: init
  type(psb_s_vect_type),intent(inout), optional   :: initu

  integer(psb_ipk_)   :: n_row,n_col
  real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act
  character           :: trans_
  character(len=20)   :: name='s_diag_solver_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select
  !
  ! For non-iterative solvers, init and initu are ignored.
  !

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()
  if (x%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,& 
         & i_err=(/itwo,n_row,izero,izero,izero/))
    goto 9999
  end if
  if (y%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,& 
         & i_err=(/ithree,n_row,izero,izero,izero/))
    goto 9999
  end if
  if (.not.allocated(sv%dv)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (sv%dv%get_nrows() < n_row) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if

  call y%mlt(alpha,sv%dv,x,beta,info,conjgx=trans_)

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='vect%mlt')
    goto 9999      
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_s_diag_solver_apply_vect

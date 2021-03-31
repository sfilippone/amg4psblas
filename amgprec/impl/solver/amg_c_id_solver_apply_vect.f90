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
subroutine amg_c_id_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
     & trans,work,wv,info,init,initu)

  use psb_base_mod
  use amg_c_id_solver, amg_protect_name => amg_c_id_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)            :: desc_data
  class(amg_c_id_solver_type), intent(inout) :: sv
  type(psb_c_vect_type),intent(inout)        :: x
  type(psb_c_vect_type),intent(inout)        :: y
  complex(psb_spk_),intent(in)                   :: alpha,beta
  character(len=1),intent(in)                  :: trans
  complex(psb_spk_),target, intent(inout)        :: work(:)
  type(psb_c_vect_type),intent(inout)        :: wv(:)
  integer(psb_ipk_), intent(out)               :: info
  character, intent(in), optional                :: init
  type(psb_c_vect_type),intent(inout), optional   :: initu

  integer(psb_ipk_)   :: n_row,n_col
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act
  character           :: trans_
  character(len=20)   :: name='c_id_solver_apply_vect'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T')
  case('C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select
  !
  ! For non-iterative solvers, init and initu are ignored.
  !

  call psb_geaxpby(alpha,x,beta,y,desc_data,info)    

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_c_id_solver_apply_vect

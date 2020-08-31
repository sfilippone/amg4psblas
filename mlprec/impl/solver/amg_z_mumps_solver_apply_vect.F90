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
subroutine z_mumps_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
     & trans,work,wv,info,init,initu)
  use psb_base_mod
  use amg_z_mumps_solver
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(amg_z_mumps_solver_type), intent(inout) :: sv
  type(psb_z_vect_type),intent(inout)  :: x
  type(psb_z_vect_type),intent(inout)  :: y
  complex(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)           :: trans
  complex(psb_dpk_),target, intent(inout) :: work(:)
  type(psb_z_vect_type),intent(inout) :: wv(:)
  integer(psb_ipk_), intent(out)       :: info
  character, intent(in), optional                 :: init
  type(psb_z_vect_type),intent(inout), optional :: initu

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='z_mumps_solver_apply_vect'

#if defined(HAVE_MUMPS_) 

  call psb_erractionsave(err_act)

  info = psb_success_
  !
  ! For non-iterative solvers, init and initu are ignored.
  !

  call x%v%sync()
  call y%v%sync()
  call sv%apply(alpha,x%v%v,beta,y%v%v,desc_data,trans,work,info)
  call y%v%set_host()
  if (info /= 0) goto 9999

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

end subroutine z_mumps_solver_apply_vect


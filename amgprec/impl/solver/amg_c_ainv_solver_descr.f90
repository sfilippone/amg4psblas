!  
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
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
subroutine amg_c_ainv_solver_descr(sv,info,iout,coarse)


  use psb_base_mod
  use amg_c_ainv_solver, amg_protect_name => amg_c_ainv_solver_descr

  Implicit None

  ! Arguments
  class(amg_c_ainv_solver_type), intent(in) :: sv
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: iout
  logical, intent(in), optional             :: coarse

  ! Local variables
  integer(psb_ipk_)      :: err_act
  integer(psb_ipk_)      :: ictxt, me, np
  character(len=20), parameter :: name='amg_c_ainv_solver_descr'
  integer(psb_ipk_)      :: iout_

  call psb_erractionsave(err_act)
  info = psb_success_
  if (present(iout)) then
    iout_ = iout
  else
    iout_ = 6
  endif

  write(iout_,*) '  AINV: Approximate Inverse with sparse biconjugation '
  write(iout_,*) '  Algorithm variant      : ',sv%algname(sv%alg)
  write(iout_,*) '  Fill level             : ',sv%fill_in
  write(iout_,*) '  Fill threshold         : ',sv%thresh

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine amg_c_ainv_solver_descr

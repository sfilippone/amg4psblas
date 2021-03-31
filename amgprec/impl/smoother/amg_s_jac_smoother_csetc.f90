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
subroutine amg_s_jac_smoother_csetc(sm,what,val,info,idx)

  use psb_base_mod
  use amg_s_jac_smoother, amg_protect_nam => amg_s_jac_smoother_csetc
  Implicit None
  ! Arguments
  class(amg_s_jac_smoother_type), intent(inout) :: sm
  character(len=*), intent(in)                     :: what
  character(len=*), intent(in)                   :: val
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), intent(in), optional       :: idx
  integer(psb_ipk_) :: err_act, ival
  character(len=20) :: name='s_jac_smoother_csetc'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(trim(what)))
    case('SMOOTHER_STOP')
      select case(psb_toupper(trim(val)))
      case('T','TRUE')
        sm%checkres   = .true.
      case('F','FALSE')
        sm%checkres   = .false.
      case default
        write(0,*) 'Unknown value for smoother_stop : "',psb_toupper(trim(val)),'"'
      end select
    case('SMOOTHER_TRACE')
      select case(psb_toupper(trim(val)))
      case('T','TRUE')
        sm%printres   = .true.
      case('F','FALSE')
        sm%printres   = .false.
      case default
        write(0,*) 'Unknown value for smoother_trace : "',psb_toupper(trim(val)),'"'
      end select
    case default
      call sm%amg_s_base_smoother_type%set(what,val,info,idx=idx)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info, name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_s_jac_smoother_csetc

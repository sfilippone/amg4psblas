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
subroutine amg_d_base_onelev_csetr(lv,what,val,info,pos,idx)
  
  use psb_base_mod
  use amg_d_onelev_mod, amg_protect_name => amg_d_base_onelev_csetr

  Implicit None

  ! Arguments
  class(amg_d_onelev_type), intent(inout) :: lv 
  character(len=*), intent(in)              :: what 
  real(psb_dpk_), intent(in)                 :: val
  integer(psb_ipk_), intent(out)            :: info
  character(len=*), optional, intent(in)    :: pos
  integer(psb_ipk_), intent(in), optional   :: idx
  ! Local 
  integer(psb_ipk_)  :: ipos_, err_act
  character(len=20) :: name='d_base_onelev_csetr'

  call psb_erractionsave(err_act)


  info = psb_success_

  select case (psb_toupper(what))

  case ('AGGR_OMEGA_VAL')
    lv%parms%aggr_omega_val= val

  case ('AGGR_THRESH')
    lv%parms%aggr_thresh   = val

  case default
    
    if (present(pos)) then
      select case(psb_toupper(trim(pos)))
      case('PRE')
        ipos_ = amg_smooth_pre_
      case('POST')
        ipos_ = amg_smooth_post_
      case default
        ipos_ = amg_smooth_both_
      end select
    else
      ipos_ = amg_smooth_both_
    end if

    if ((ipos_==amg_smooth_pre_) .or.(ipos_==amg_smooth_both_)) then 
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info,idx=idx)
      end if
    end if
    if ((ipos_==amg_smooth_post_).or.(ipos_==amg_smooth_both_))then 
      if (allocated(lv%sm2a)) then 
        call lv%sm2a%set(what,val,info,idx=idx)
      end if
    end if
    if (allocated(lv%aggr)) call lv%aggr%set(what,val,info,idx=idx)

  end select

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_base_onelev_csetr

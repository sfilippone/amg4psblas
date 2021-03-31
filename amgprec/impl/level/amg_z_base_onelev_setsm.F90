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
subroutine amg_z_base_onelev_setsm(lev,val,info,pos)

  use psb_base_mod
  use amg_z_prec_mod, amg_protect_name => amg_z_base_onelev_setsm

  implicit none

  ! Arguments
  class(amg_z_onelev_type), target, intent(inout) :: lev
  class(amg_z_base_smoother_type), intent(in)  :: val
  integer(psb_ipk_), intent(out)               :: info
  character(len=*), optional, intent(in)       :: pos
  
  ! Local variables
  integer(psb_ipk_)                      :: ipos_
  character(len=*), parameter            :: name='amg_base_onelev_setsm'
  
  info = psb_success_

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

  if (ipos_ == amg_smooth_both_) then
    if (allocated(lev%sm2a)) then
      call lev%sm2a%free(info)
      deallocate(lev%sm2a, stat=info)
      lev%sm2 => null()
    end if
  end if
     
  select case(ipos_)
  case(amg_smooth_pre_, amg_smooth_both_) 
    if (allocated(lev%sm)) then
      if (.not.same_type_as(lev%sm,val)) then
        call lev%sm%free(info)
        deallocate(lev%sm, stat=info)
      end if
    endif
    if (.not.allocated(lev%sm)) then
#ifdef HAVE_MOLD 
      allocate(lev%sm,mold=val) 
#else
      allocate(lev%sm,source=val) 
#endif
    end if
    call lev%sm%default()        
    if (ipos_ ==  amg_smooth_both_) lev%sm2 => lev%sm
  case(amg_smooth_post_) 
    if (allocated(lev%sm2a)) then
      if (.not.same_type_as(lev%sm2a,val)) then
        call lev%sm2a%free(info)
        deallocate(lev%sm2a, stat=info)
      endif
    end if
    if (.not.allocated(lev%sm2a)) then
#ifdef HAVE_MOLD 
      allocate(lev%sm2a,mold=val) 
#else
      allocate(lev%sm2a,source=val) 
#endif
    end if
    call lev%sm2a%default()
    lev%sm2 => lev%sm2a
  end select
    
end subroutine amg_z_base_onelev_setsm


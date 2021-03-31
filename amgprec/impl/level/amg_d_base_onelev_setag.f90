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
subroutine amg_d_base_onelev_setag(lv,val,info,pos)

  use psb_base_mod
  use amg_d_onelev_mod, amg_protect_name => amg_d_base_onelev_setag

  implicit none

  ! Arguments
  class(amg_d_onelev_type), target, intent(inout) :: lv
  class(amg_d_base_aggregator_type), intent(in)   :: val
  integer(psb_ipk_), intent(out)                  :: info
  character(len=*), optional, intent(in)          :: pos
  
  ! Local variables
  integer(psb_ipk_)                :: ipos_
  character(len=*), parameter      :: name='amg_base_onelev_setag'

  info = psb_success_

  ! Ignore pos for aggregator
  
  if (allocated(lv%aggr)) then 
    if (.not.same_type_as(lv%aggr,val))  then
      call lv%aggr%free(info)
      deallocate(lv%aggr,stat=info)
      if (info /= 0) then
        info = 3111
        return
      end if
    end if
  end if
      
  if (.not.allocated(lv%aggr)) then 
    allocate(lv%aggr,mold=val,stat=info) 
    if (info /= 0) then
      info = 3111
      return
    end if
    lv%parms%par_aggr_alg  = amg_ext_aggr_
    lv%parms%aggr_type     = amg_noalg_
    call lv%aggr%default()
  end if
  
end subroutine amg_d_base_onelev_setag


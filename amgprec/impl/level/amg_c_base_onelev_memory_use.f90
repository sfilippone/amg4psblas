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
!
!
! verbosity:
!        <0: suppress all messages
!         0: normal
!        >1: increased details 
!
subroutine amg_c_base_onelev_memory_use(lv,il,nl,ilmin,info,iout,verbosity,prefix)
  
  use psb_base_mod
  use amg_c_onelev_mod, amg_protect_name => amg_c_base_onelev_memory_use
  Implicit None
  ! Arguments
  class(amg_c_onelev_type), intent(in)    :: lv
  integer(psb_ipk_), intent(in)           :: il,nl,ilmin
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), intent(in), optional :: iout
  integer(psb_ipk_), intent(in), optional :: verbosity
  character(len=*), intent(in), optional  :: prefix


  ! Local variables
  integer(psb_ipk_)  :: err_act
  character(len=20), parameter :: name='amg_c_base_onelev_memory_use'
  integer(psb_ipk_)  :: iout_, verbosity_
  logical      :: coarse
  character(1024)    :: prefix_


  call psb_erractionsave(err_act)


  coarse = (il==nl)

  if (present(iout)) then 
    iout_ = iout
  else 
    iout_ = psb_out_unit
  end if
  
  if (present(verbosity)) then
    verbosity_ = verbosity
  else
    verbosity_ = 0
  end if
  if (verbosity_ < 0) goto 9998
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if

  write(iout_,*) trim(prefix_)
  
  if (coarse)  then 
    write(iout_,*) trim(prefix_), ' Level ',il,' (coarse)'
  else
    write(iout_,*) trim(prefix_), ' Level ',il
  end if
  
  write(iout_,*) trim(prefix_), '           Matrix:', lv%base_a%sizeof()
  write(iout_,*) trim(prefix_), '       Descriptor:', lv%base_desc%sizeof()
  if (il >1) write(iout_,*) trim(prefix_), '       Linear map:', lv%linmap%sizeof()
  if (allocated(lv%sm)) write(iout_,*) trim(prefix_), '         Smoother:', lv%sm%sizeof()
  if (allocated(lv%sm2a)) write(iout_,*) trim(prefix_), '       Smoother 2:', lv%sm2a%sizeof()
  if (allocated(lv%wrk)) write(iout_,*) trim(prefix_), '        Workspace:', lv%wrk%sizeof()
  

9998 continue
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_c_base_onelev_memory_use

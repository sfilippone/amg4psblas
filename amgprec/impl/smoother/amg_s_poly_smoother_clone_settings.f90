!  
!   
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               asd on PSBLAS (Parallel Sparse BLAS version 3.7)
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
subroutine amg_s_poly_smoother_clone_settings(sm,smout,info)

  use psb_base_mod
  use amg_s_poly_smoother, amg_protect_name =>  amg_s_poly_smoother_clone_settings
  Implicit None
  ! Arguments
  class(amg_s_poly_smoother_type), intent(inout)              :: sm
  class(amg_s_base_smoother_type), intent(inout) :: smout
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_)  :: err_act
  character(len=20) :: name='d_poly_smoother_clone_settings'

  call psb_erractionsave(err_act)

  info = psb_success_
  
  select type(smout)
  class is(amg_s_poly_smoother_type)

    smout%pa => null()
    smout%pdegree   = sm%pdegree
    smout%variant   = sm%variant
    smout%cf_a      = sm%cf_a
    smout%rho_ba    = sm%rho_ba
    smout%rho_estimate  = sm%rho_estimate
    smout%rho_estimate_iterations  = sm%rho_estimate_iterations
    smout%poly_beta = sm%poly_beta
    
    if (allocated(smout%sv)) then
      if (.not.same_type_as(sm%sv,smout%sv)) then
        call smout%sv%free(info)
        if (info == 0) deallocate(smout%sv,stat=info)
      end if
    end if
    if (info /= 0) then
      info = psb_err_internal_error_
    else
      if (allocated(smout%sv)) then
        if (same_type_as(sm%sv,smout%sv)) then
          call sm%sv%clone_settings(smout%sv,info)
        else
          info = psb_err_internal_error_
        end if
      else
        allocate(smout%sv,mold=sm%sv,stat=info)
        if (info == 0) call sm%sv%clone_settings(smout%sv,info)
        if (info /= 0) info = psb_err_internal_error_
      end if
    end if
  class default
    info = psb_err_internal_error_
  end select
  
  if (info /= 0) then
    call psb_errpush(info,name) 
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_s_poly_smoother_clone_settings

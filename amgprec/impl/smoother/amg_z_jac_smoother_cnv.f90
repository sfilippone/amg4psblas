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
subroutine amg_z_jac_smoother_cnv(sm,info,amold,vmold,imold)

  use psb_base_mod
  use amg_z_diag_solver
  use amg_z_jac_smoother, amg_protect_name => amg_z_jac_smoother_cnv
  Implicit None

  ! Arguments
  class(amg_z_jac_smoother_type), intent(inout)      :: sm
  integer(psb_ipk_), intent(out)                     :: info
  class(psb_z_base_sparse_mat), intent(in), optional :: amold
  class(psb_z_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  ! Local variables
  integer(psb_ipk_) :: i, err_act, debug_unit, debug_level
  character(len=20) :: name='d_jac_smoother_cnv', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  if (info == psb_success_) then
    if (sm%nd%is_asb()) then
      if (present(amold)) then 
        call sm%nd%cscnv(info,&
             & mold=amold,dupl=psb_dupl_add_)
      else
        call sm%nd%cscnv(info,&
             & type='csr',dupl=psb_dupl_add_)
      endif
    end if
  end if

  if  (info == psb_success_) then 
    if (allocated(sm%sv)) &
         & call sm%sv%cnv(info,amold=amold,vmold=vmold,imold=imold)
  end if
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='solver cnv')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_z_jac_smoother_cnv

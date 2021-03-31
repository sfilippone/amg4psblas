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
subroutine amg_s_base_onelev_cnv(lv,info,amold,vmold,imold)
  
  use psb_base_mod
  use amg_s_onelev_mod, amg_protect_name =>  amg_s_base_onelev_cnv
  implicit none 

  class(amg_s_onelev_type), intent(inout) :: lv
  integer(psb_ipk_), intent(out)          :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  integer(psb_ipk_) :: i

  info = psb_success_
  
  if (any((/present(amold),present(vmold),present(imold)/))) then 
    if (allocated(lv%sm)) &
         & call lv%sm%cnv(info,amold=amold,vmold=vmold,imold=imold)
    if (info == psb_success_ .and. allocated(lv%sm2a)) &
         & call lv%sm2a%cnv(info,amold=amold,vmold=vmold,imold=imold)
    if (info == psb_success_ .and. allocated(lv%wrk)) &
         & call lv%wrk%cnv(info,vmold=vmold)
    if (info == psb_success_.and. lv%ac%is_asb()) &
         & call lv%ac%cscnv(info,mold=amold)
    if (info == psb_success_ .and. lv%desc_ac%is_ok() &
         & .and. present(imold)) call lv%desc_ac%cnv(imold)
    if (info == psb_success_) call lv%linmap%cnv(info,mold=amold,imold=imold)
  end if
end subroutine amg_s_base_onelev_cnv

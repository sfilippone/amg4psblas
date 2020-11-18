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
subroutine amg_s_as_smoother_cnv(sm,info,amold,vmold,imold)
  
  use psb_base_mod
  use amg_s_as_smoother, amg_protect_nam => amg_s_as_smoother_cnv
  Implicit None

  ! Arguments
  class(amg_s_as_smoother_type), intent(inout)     :: sm
  integer(psb_ipk_), intent(out)                     :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  ! Local variables
  type(psb_dspmat_type) :: blck, atmp
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nhalo, novr, data_, nzeros
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act, debug_unit, debug_level
  character(len=20)   :: name='d_as_smoother_cnv', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = sm%desc_data%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  if (allocated(sm%sv)) &
       & call sm%sv%cnv(info,amold=amold,vmold=vmold,imold=imold)

  if (info == psb_success_) then 
    if (present(amold)) then 
      if (sm%nd%is_asb()) call sm%nd%cscnv(info,&
           & mold=amold,dupl=psb_dupl_add_)
    end if
  end if
  if (info == psb_success_) then
    if (present(imold)) then 
      call sm%desc_data%cnv(imold)
    end if
  end if
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_s_as_smoother_cnv

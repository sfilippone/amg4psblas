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
subroutine amg_d_ainv_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)


  use psb_base_mod
  use psb_d_ainv_fact_mod
  use amg_d_ainv_solver, amg_protect_name => amg_d_ainv_solver_bld

  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target    :: a
  Type(psb_desc_type), Intent(inout)           :: desc_a
  class(amg_d_ainv_solver_type), intent(inout) :: sv
  integer(psb_ipk_), intent(out)                      :: info
  type(psb_dspmat_type), intent(in), target, optional :: b
  class(psb_d_base_sparse_mat), intent(in), optional  :: amold
  class(psb_d_base_vect_type), intent(in), optional   :: vmold
  class(psb_i_base_vect_type), intent(in), optional   :: imold
  ! Local variables
  integer(psb_ipk_)  :: n_row,n_col, nrow_a, nztota
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer(psb_ipk_)  :: np,me,i, err_act, debug_unit, debug_level
  type(psb_ctxt_type) :: ctxt
  character(len=20)  :: name='amg_d_ainv_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  call psb_d_ainv_bld(a,sv%alg,sv%fill_in,sv%thresh,&
       & sv%w,sv%d,sv%z,desc_a,info,b,iscale=amg_ilu_scale_maxval_)

  if ((info == psb_success_) .and.present(amold)) then
    call sv%w%set_asb()
    call sv%w%trim()
    call sv%z%set_asb()
    call sv%z%trim()
    call sv%w%cscnv(info,mold=amold)
    if (info == psb_success_) &
         & call sv%z%cscnv(info,mold=amold)
  end if

  if (info == psb_success_) &
       & call sv%dv%bld(sv%d,mold=vmold)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name)
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine amg_d_ainv_solver_bld

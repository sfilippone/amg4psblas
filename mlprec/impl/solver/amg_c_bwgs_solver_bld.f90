!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
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
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine amg_c_bwgs_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

  use psb_base_mod
  use amg_c_gs_solver, amg_protect_name => amg_c_bwgs_solver_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(in)                     :: desc_a 
  class(amg_c_bwgs_solver_type), intent(inout)         :: sv
  integer(psb_ipk_), intent(out)                      :: info
  type(psb_cspmat_type), intent(in), target, optional :: b
  class(psb_c_base_sparse_mat), intent(in), optional  :: amold
  class(psb_c_base_vect_type), intent(in), optional   :: vmold
  class(psb_i_base_vect_type), intent(in), optional   :: imold
  ! Local variables
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota
  integer(psb_ipk_) :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20) :: name='d_bwgs_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()

  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()

  if (sv%eps <= dzero) then
    !
    ! This cuts out the off-diagonal part, because it's supposed to
    ! be handled by the outer Jacobi smoother.
    !
    call a%tril(sv%l,info,diag=-ione,jmax=nrow_a,u=sv%u)     

  else

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999       
  end if



  call sv%l%set_asb()
  call sv%l%trim()
  call sv%u%set_asb()
  call sv%u%trim()

  if (present(amold)) then 
    call sv%l%cscnv(info,mold=amold)
    call sv%u%cscnv(info,mold=amold)
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_c_bwgs_solver_bld

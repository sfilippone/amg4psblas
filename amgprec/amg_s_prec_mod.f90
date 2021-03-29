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
! File: amg_s_prec_mod.f90
!
! Module: amg_s_prec_mod
!
!  This module defines the user interfaces to the real/complex, single/double
!  precision versions of the user-level MLD2P4 routines.
!
module amg_s_prec_mod

  use amg_s_prec_type
  use amg_s_jac_smoother
  use amg_s_as_smoother
  use amg_s_id_solver
  use amg_s_diag_solver
  use amg_s_l1_diag_solver
  use amg_s_ilu_solver
  use amg_s_gs_solver
  use amg_s_ainv_solver
  use amg_s_invk_solver
  use amg_s_invt_solver
  use amg_s_krm_solver

  interface amg_extprol_bld
    subroutine amg_s_extprol_bld(a,desc_a,p,prolv,restrv,info,amold,vmold,imold)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & psb_s_base_sparse_mat, psb_s_base_vect_type, &
           & psb_i_base_vect_type, amg_sprec_type, psb_ipk_

      ! Arguments
      type(psb_sspmat_type),intent(in), target           :: a
      type(psb_sspmat_type),intent(inout), target        :: prolv(:)
      type(psb_sspmat_type),intent(inout), target        :: restrv(:)
      type(psb_desc_type), intent(inout), target         :: desc_a
      type(amg_sprec_type),intent(inout),target          :: p
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: amold
      class(psb_s_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      ! !$  character, intent(in), optional         :: upd
    end subroutine amg_s_extprol_bld
  end interface amg_extprol_bld

end module amg_s_prec_mod

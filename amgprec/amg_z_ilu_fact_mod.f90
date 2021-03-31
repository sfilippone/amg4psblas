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
! File: amg_z_ilu_fact_mod.f90
!
! Module: amg_z_ilu_fact_mod
!
!  This module defines some interfaces used internally by the implementation of
!  amg_z_ilu_solver, but not visible to the end user. 
!
!
module amg_z_ilu_fact_mod

  use amg_z_base_solver_mod 

  interface amg_ilu0_fact
    subroutine amg_zilu0_fact(ialg,a,l,u,d,info,blck,upd)
      import psb_zspmat_type, psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)         :: ialg
      integer(psb_ipk_), intent(out)        :: info
      type(psb_zspmat_type),intent(in)      :: a
      type(psb_zspmat_type),intent(inout)   :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      character, intent(in), optional       :: upd
      complex(psb_dpk_), intent(inout)      :: d(:)
    end subroutine amg_zilu0_fact
  end interface

  interface amg_iluk_fact
    subroutine amg_ziluk_fact(fill_in,ialg,a,l,u,d,info,blck)
      import psb_zspmat_type, psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)        :: fill_in,ialg
      integer(psb_ipk_), intent(out)       :: info
      type(psb_zspmat_type),intent(in)     :: a
      type(psb_zspmat_type),intent(inout)  :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine amg_ziluk_fact
  end interface

  interface amg_ilut_fact
    subroutine amg_zilut_fact(fill_in,thres,a,l,u,d,info,blck,iscale)
      import  psb_zspmat_type, psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)        :: fill_in
      real(psb_dpk_), intent(in)           :: thres
      integer(psb_ipk_), intent(out)       :: info
      type(psb_zspmat_type),intent(in)     :: a
      type(psb_zspmat_type),intent(inout)  :: l,u
      complex(psb_dpk_), intent(inout)     :: d(:)
      type(psb_zspmat_type),intent(in), optional, target :: blck
      integer(psb_ipk_), intent(in), optional  :: iscale
    end subroutine amg_zilut_fact
  end interface

end module amg_z_ilu_fact_mod

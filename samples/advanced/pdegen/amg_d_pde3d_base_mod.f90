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
module amg_d_pde3d_base_mod
  use psb_base_mod, only : psb_dpk_, done, dzero
  real(psb_dpk_), save, private :: epsilon=done/80
contains
  subroutine pde_set_parm3d_base(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm3d_base
  !
  ! functions parametrizing the differential equation
  !
  function b1_base(x,y,z)
    implicit none 
    real(psb_dpk_) :: b1_base
    real(psb_dpk_), intent(in) :: x,y,z
    b1_base=dzero/sqrt(3.0_psb_dpk_)
  end function b1_base
  function b2_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  b2_base
    real(psb_dpk_), intent(in) :: x,y,z
    b2_base=dzero/sqrt(3.0_psb_dpk_)
  end function b2_base
  function b3_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  b3_base
    real(psb_dpk_), intent(in) :: x,y,z
    b3_base=dzero/sqrt(3.0_psb_dpk_)
  end function b3_base
  function c_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  c_base
    real(psb_dpk_), intent(in) :: x,y,z
    c_base=dzero
  end function c_base
  function a1_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a1_base
    real(psb_dpk_), intent(in) :: x,y,z
    a1_base=epsilon
  end function a1_base
  function a2_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a2_base
    real(psb_dpk_), intent(in) :: x,y,z
    a2_base=epsilon
  end function a2_base
  function a3_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a3_base
    real(psb_dpk_), intent(in) :: x,y,z
    a3_base=epsilon
  end function a3_base
  function g_base(x,y,z)
    implicit none 
    real(psb_dpk_) ::  g_base
    real(psb_dpk_), intent(in) :: x,y,z
    g_base = dzero
    if (x == done) then
      g_base = done
    else if (x == dzero) then
      g_base = done
    end if
  end function g_base
end module amg_d_pde3d_base_mod

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
module amg_s_pde3d_gauss_mod
  use psb_base_mod, only : psb_spk_, sone, szero
  real(psb_spk_), save, private :: epsilon=sone/80
contains
  subroutine pde_set_parm3d_gauss(dat)
    real(psb_spk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm3d_gauss
  !
  ! functions parametrizing the differential equation
  !
  function b1_gauss(x,y,z)
    implicit none 
    real(psb_spk_) :: b1_gauss
    real(psb_spk_), intent(in) :: x,y,z
    b1_gauss=sone/sqrt(3.0_psb_spk_)-2*x*exp(-(x**2+y**2+z**2))
  end function b1_gauss
  function b2_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  b2_gauss
    real(psb_spk_), intent(in) :: x,y,z
    b2_gauss=sone/sqrt(3.0_psb_spk_)-2*y*exp(-(x**2+y**2+z**2))
  end function b2_gauss
  function b3_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  b3_gauss
    real(psb_spk_), intent(in) :: x,y,z
    b3_gauss=sone/sqrt(3.0_psb_spk_)-2*z*exp(-(x**2+y**2+z**2))
  end function b3_gauss
  function c_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  c_gauss
    real(psb_spk_), intent(in) :: x,y,z
    c_gauss=szero
  end function c_gauss
  function a1_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  a1_gauss
    real(psb_spk_), intent(in) :: x,y,z
    a1_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a1_gauss
  function a2_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  a2_gauss
    real(psb_spk_), intent(in) :: x,y,z
    a2_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a2_gauss
  function a3_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  a3_gauss
    real(psb_spk_), intent(in) :: x,y,z
    a3_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a3_gauss
  function g_gauss(x,y,z)
    implicit none 
    real(psb_spk_) ::  g_gauss
    real(psb_spk_), intent(in) :: x,y,z
    g_gauss = szero
    if (x == sone) then
      g_gauss = sone
    else if (x == szero) then
      g_gauss = sone
    end if
  end function g_gauss
end module amg_s_pde3d_gauss_mod

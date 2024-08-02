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
module amg_d_pde3d_const_mod
  use psb_base_mod, only : psb_dpk_, done, dzero
  real(psb_dpk_), save, private :: epsilon=done/80
contains
  subroutine pde_set_parm3d_const(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm3d_const
  !
  ! functions parametrizing the differential equation
  !
  function b1_const(x,y,z)
    implicit none 
    real(psb_dpk_) :: b1_const
    real(psb_dpk_), intent(in) :: x,y,z
    b1_const=done/sqrt(3.0_psb_dpk_)
  end function b1_const
  function b2_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  b2_const
    real(psb_dpk_), intent(in) :: x,y,z
    b2_const=done/sqrt(3.0_psb_dpk_)
  end function b2_const
  function b3_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  b3_const
    real(psb_dpk_), intent(in) :: x,y,z
    b3_const=dzero/sqrt(3.0_psb_dpk_)
  end function b3_const
  function c_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  c_const
    real(psb_dpk_), intent(in) :: x,y,z
    c_const=dzero
  end function c_const
  function a1_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a1_const
    real(psb_dpk_), intent(in) :: x,y,z
    a1_const=epsilon
  end function a1_const
  function a2_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a2_const
    real(psb_dpk_), intent(in) :: x,y,z
    a2_const=epsilon
  end function a2_const
  function a3_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a3_const
    real(psb_dpk_), intent(in) :: x,y,z
    a3_const=epsilon
  end function a3_const
  function g_const(x,y,z)
    implicit none 
    real(psb_dpk_) ::  g_const
    real(psb_dpk_), intent(in) :: x,y,z
    g_const = dzero
    if (x == done) then
      g_const = done
    else if (x == dzero) then
      g_const = done
    end if
  end function g_const
end module amg_d_pde3d_const_mod

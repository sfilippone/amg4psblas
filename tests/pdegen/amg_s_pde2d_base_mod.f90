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
module amg_s_pde2d_base_mod
  use psb_base_mod, only : psb_spk_, szero, sone
  real(psb_spk_), save, private :: epsilon=sone/80
contains
  subroutine pde_set_parm(dat)
    real(psb_spk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y
    b1 = szero/1.414_psb_spk_
  end function b1
  function b2(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y
    b2 = szero/1.414_psb_spk_
  end function b2
  function c(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  c
    real(psb_spk_), intent(in) :: x,y
    c = szero
  end function c
  function a1(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  a1
    real(psb_spk_), intent(in) :: x,y
    a1=sone*epsilon
  end function a1
  function a2(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  a2
    real(psb_spk_), intent(in) :: x,y
    a2=sone*epsilon
  end function a2
  function g(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  g
    real(psb_spk_), intent(in) :: x,y
    g = szero
    if (x == sone) then
      g = sone
    else if (x == szero) then
      g = sone
    end if
  end function g
end module amg_s_pde2d_base_mod

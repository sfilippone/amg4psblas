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
!         software without specific prior written permission.
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
module amg_d_pde3d_aniso_mod

  use psb_base_mod, only : psb_dpk_, done, dzero

  real(psb_dpk_), save, private :: epsilon = done/80
  real(psb_dpk_), save, private :: theta   = dzero

contains

  subroutine pde_set_parm3d_aniso(dat)

    real(psb_dpk_), intent(in) :: dat(:)

    epsilon = dat(1)
    theta   = dat(2)

  end subroutine pde_set_parm3d_aniso

  !

  ! Diffusion tensor coefficients
  !
  !        [ k11  k12  k13 ]
  !  K  =  [ k12  k22  k23 ]
  !        [ k13  k23  k33 ]
  !
  ! corresponding to a rotation by theta in the x-y plane
  ! of the diagonal tensor diag(epsilon,1,1).

  !

  function k11(x,y,z)

    implicit none

    real(psb_dpk_) :: k11
    real(psb_dpk_), intent(in) :: x,y,z

    k11 = epsilon*cos(theta)**2 + sin(theta)**2

  end function k11

  function k22(x,y,z)

    implicit none

    real(psb_dpk_) :: k22
    real(psb_dpk_), intent(in) :: x,y,z

    k22 = epsilon*sin(theta)**2 + cos(theta)**2

  end function k22

  function k33(x,y,z)

    implicit none

    real(psb_dpk_) :: k33
    real(psb_dpk_), intent(in) :: x,y,z

    k33 = done

  end function k33

  function k12(x,y,z)

    implicit none

    real(psb_dpk_) :: k12
    real(psb_dpk_), intent(in) :: x,y,z

    k12 = (epsilon-done)*sin(theta)*cos(theta)

  end function k12

  function k13(x,y,z)

    implicit none

    real(psb_dpk_) :: k13
    real(psb_dpk_), intent(in) :: x,y,z

    k13 = dzero

  end function k13

  function k23(x,y,z)

    implicit none

    real(psb_dpk_) :: k23
    real(psb_dpk_), intent(in) :: x,y,z

    k23 = dzero

  end function k23

  !

  ! Right-hand side / boundary data

  !

  function g_aniso(x,y,z)

    implicit none

    real(psb_dpk_) :: g_aniso
    real(psb_dpk_), intent(in) :: x,y,z

    g_aniso = dzero

    if (x == done) then

      g_aniso = done

    else if (x == dzero) then

      g_aniso = done

    end if

  end function g_aniso

end module amg_d_pde3d_aniso_mod
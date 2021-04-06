module amg_d_pde3d_gauss_mod
  use psb_base_mod, only : psb_dpk_, done
  real(psb_dpk_), save, private :: epsilon=done/80
contains
  subroutine pde_set_parm(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) :: b1_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    b1_gauss=done/sqrt(3.0_psb_dpk_)-2*x*exp(-(x**2+y**2+z**2))
  end function b1_gauss
  function b2_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  b2_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    b2_gauss=done/sqrt(3.0_psb_dpk_)-2*y*exp(-(x**2+y**2+z**2))
  end function b2_gauss
  function b3_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  b3_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    b3_gauss=done/sqrt(3.0_psb_dpk_)-2*z*exp(-(x**2+y**2+z**2))
  end function b3_gauss
  function c_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_, dzero
    real(psb_dpk_) ::  c_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    c=dzero
  end function c_gauss
  function a1_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a1_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    a1_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a1_gauss
  function a2_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a2_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    a2_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a2_gauss
  function a3_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a3_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    a3_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a3_gauss
  function g_gauss(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  g_gauss
    real(psb_dpk_), intent(in) :: x,y,z
    g_gauss = dzero
    if (x == done) then
      g_gauss = done
    else if (x == dzero) then
      g_gauss = done
    end if
  end function g_gauss
end module amg_d_pde3d_gauss_mod

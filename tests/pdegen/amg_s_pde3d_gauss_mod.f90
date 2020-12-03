module amg_s_pde3d_gauss_mod
  use psb_base_mod, only : psb_spk_, sone
  real(psb_spk_), save, private :: epsilon=sone/80
contains
  subroutine pde_set_parm(dat)
    real(psb_spk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) :: b1_gauss
    real(psb_spk_), intent(in) :: x,y,z
    b1_gauss=sone/sqrt(3.0_psb_spk_)-2*x*exp(-(x**2+y**2+z**2))
  end function b1_gauss
  function b2_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) ::  b2_gauss
    real(psb_spk_), intent(in) :: x,y,z
    b2_gauss=sone/sqrt(3.0_psb_spk_)-2*y*exp(-(x**2+y**2+z**2))
  end function b2_gauss
  function b3_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) ::  b3_gauss
    real(psb_spk_), intent(in) :: x,y,z
    b3_gauss=sone/sqrt(3.0_psb_spk_)-2*z*exp(-(x**2+y**2+z**2))
  end function b3_gauss
  function c_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_, szero
    real(psb_spk_) ::  c_gauss
    real(psb_spk_), intent(in) :: x,y,z
    c=szero
  end function c_gauss
  function a1_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a1_gauss
    real(psb_spk_), intent(in) :: x,y,z
    a1_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a1_gauss
  function a2_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a2_gauss
    real(psb_spk_), intent(in) :: x,y,z
    a2_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a2_gauss
  function a3_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a3_gauss
    real(psb_spk_), intent(in) :: x,y,z
    a3_gauss=epsilon*exp(-(x**2+y**2+z**2))
  end function a3_gauss
  function g_gauss(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
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

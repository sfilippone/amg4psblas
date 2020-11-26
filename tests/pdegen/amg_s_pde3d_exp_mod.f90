module amg_s_pde3d_exp_mod
  use psb_base_mod, only : psb_spk_, sone
  real(psb_spk_), save, private :: epsilon=sone/160
contains
  subroutine pde_set_parm(dat)
    real(psb_spk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_exp(x,y,z)
    use psb_base_mod, only : psb_spk_, szero
    real(psb_spk_) :: b1_exp
    real(psb_spk_), intent(in) :: x,y,z
    b1_exp=szero/sqrt(3.0_psb_spk_)
  end function b1_exp
  function b2_exp(x,y,z)
    use psb_base_mod, only : psb_spk_, szero
    real(psb_spk_) ::  b2_exp
    real(psb_spk_), intent(in) :: x,y,z
    b2_exp=szero/sqrt(3.0_psb_spk_)
  end function b2_exp
  function b3_exp(x,y,z)
    use psb_base_mod, only : psb_spk_, szero
    real(psb_spk_) ::  b3_exp
    real(psb_spk_), intent(in) :: x,y,z
    b3_exp=szero/sqrt(3.0_psb_spk_)
  end function b3_exp
  function c_exp(x,y,z)
    use psb_base_mod, only : psb_spk_, szero
    real(psb_spk_) ::  c_exp
    real(psb_spk_), intent(in) :: x,y,z
    c_exp=szero
  end function c_exp
  function a1_exp(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a1_exp
    real(psb_spk_), intent(in) :: x,y,z
    a1_exp=epsilon*exp(-(x+y+z))
  end function a1_exp
  function a2_exp(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a2_exp
    real(psb_spk_), intent(in) :: x,y,z
    a2_exp=epsilon*exp(-(x+y+z))
  end function a2_exp
  function a3_exp(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a3_exp
    real(psb_spk_), intent(in) :: x,y,z
    a3_exp=epsilon*exp(-(x+y+z))
  end function a3_exp
  function g_exp(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  g_exp
    real(psb_spk_), intent(in) :: x,y,z
    g_exp = szero
    if (x == sone) then
      g_exp = sone
    else if (x == szero) then
      g_exp = sone
    end if
  end function g_exp
end module amg_s_pde3d_exp_mod

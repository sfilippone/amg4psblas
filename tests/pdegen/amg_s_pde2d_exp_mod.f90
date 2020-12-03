module amg_s_pde2d_exp_mod
  use psb_base_mod, only : psb_spk_, sone, szero
  real(psb_spk_), save, private :: epsilon=sone/80
contains
  subroutine pde_set_parm(dat)
    real(psb_spk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_exp(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) :: b1_exp
    real(psb_spk_), intent(in) :: x,y
    b1_exp = szero
  end function b1_exp
  function b2_exp(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  b2_exp
    real(psb_spk_), intent(in) :: x,y
    b2_exp = szero
  end function b2_exp
  function c_exp(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  c_exp
    real(psb_spk_), intent(in) :: x,y
    c_exp = szero
  end function c_exp
  function a1_exp(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  a1_exp
    real(psb_spk_), intent(in) :: x,y
    a1=sone*epsilon*exp(-(x+y))
  end function a1_exp
  function a2_exp(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  a2_exp
    real(psb_spk_), intent(in) :: x,y
    a2=sone*epsilon*exp(-(x+y))
  end function a2_exp
  function g_exp(x,y)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  g_exp
    real(psb_spk_), intent(in) :: x,y
    g_exp = szero
    if (x == sone) then
      g_exp = sone
    else if (x == szero) then
      g_exp = sone
    end if
  end function g_exp
end module amg_s_pde2d_exp_mod

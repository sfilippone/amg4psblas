module amg_d_pde2d_exp_mod
  use psb_base_mod, only : psb_dpk_, done, dzero
  real(psb_dpk_), save, private :: epsilon=done/80
contains
  subroutine pde_set_parm(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_exp(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) :: b1_exp
    real(psb_dpk_), intent(in) :: x,y
    b1_exp = dzero
  end function b1_exp
  function b2_exp(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  b2_exp
    real(psb_dpk_), intent(in) :: x,y
    b2_exp = dzero
  end function b2_exp
  function c_exp(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  c_exp
    real(psb_dpk_), intent(in) :: x,y
    c_exp = dzero
  end function c_exp
  function a1_exp(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  a1_exp
    real(psb_dpk_), intent(in) :: x,y
    a1=done*epsilon*exp(-(x+y))
  end function a1_exp
  function a2_exp(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  a2_exp
    real(psb_dpk_), intent(in) :: x,y
    a2=done*epsilon*exp(-(x+y))
  end function a2_exp
  function g_exp(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  g_exp
    real(psb_dpk_), intent(in) :: x,y
    g_exp = dzero
    if (x == done) then
      g_exp = done
    else if (x == dzero) then
      g_exp = done
    end if
  end function g_exp
end module amg_d_pde2d_exp_mod

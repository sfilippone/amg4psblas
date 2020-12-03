module amg_d_pde3d_exp_mod
  use psb_base_mod, only : psb_dpk_, done
  real(psb_dpk_), save, private :: epsilon=done/160
contains
  subroutine pde_set_parm(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_, dzero
    real(psb_dpk_) :: b1_exp
    real(psb_dpk_), intent(in) :: x,y,z
    b1_exp=dzero/sqrt(3.0_psb_dpk_)
  end function b1_exp
  function b2_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_, dzero
    real(psb_dpk_) ::  b2_exp
    real(psb_dpk_), intent(in) :: x,y,z
    b2_exp=dzero/sqrt(3.0_psb_dpk_)
  end function b2_exp
  function b3_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_, dzero
    real(psb_dpk_) ::  b3_exp
    real(psb_dpk_), intent(in) :: x,y,z
    b3_exp=dzero/sqrt(3.0_psb_dpk_)
  end function b3_exp
  function c_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_, dzero
    real(psb_dpk_) ::  c_exp
    real(psb_dpk_), intent(in) :: x,y,z
    c_exp=dzero
  end function c_exp
  function a1_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a1_exp
    real(psb_dpk_), intent(in) :: x,y,z
    a1_exp=epsilon*exp(-(x+y+z))
  end function a1_exp
  function a2_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a2_exp
    real(psb_dpk_), intent(in) :: x,y,z
    a2_exp=epsilon*exp(-(x+y+z))
  end function a2_exp
  function a3_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a3_exp
    real(psb_dpk_), intent(in) :: x,y,z
    a3_exp=epsilon*exp(-(x+y+z))
  end function a3_exp
  function g_exp(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  g_exp
    real(psb_dpk_), intent(in) :: x,y,z
    g_exp = dzero
    if (x == done) then
      g_exp = done
    else if (x == dzero) then
      g_exp = done
    end if
  end function g_exp
end module amg_d_pde3d_exp_mod

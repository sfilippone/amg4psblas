module amg_d_pde2d_box_mod
  use psb_base_mod, only : psb_dpk_, dzero, done
  real(psb_dpk_), save, private :: epsilon=done/80
contains
  subroutine pde_set_parm(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1_box(x,y)
    use psb_base_mod, only : psb_dpk_, dzero, done
    real(psb_dpk_) :: b1_box
    real(psb_dpk_), intent(in) :: x,y
    b1_box = done/1.414_psb_dpk_
  end function b1_box
  function b2_box(x,y)
    use psb_base_mod, only : psb_dpk_, dzero, done
    real(psb_dpk_) ::  b2_box
    real(psb_dpk_), intent(in) :: x,y
    b2_box = done/1.414_psb_dpk_
  end function b2_box
  function c_box(x,y)
    use psb_base_mod, only : psb_dpk_, dzero, done
    real(psb_dpk_) ::  c_box
    real(psb_dpk_), intent(in) :: x,y
    c_box = dzero
  end function c_box
  function a1_box(x,y)
    use psb_base_mod, only : psb_dpk_, dzero, done
    real(psb_dpk_) ::  a1_box
    real(psb_dpk_), intent(in) :: x,y
    a1_box=done*epsilon
  end function a1_box
  function a2_box(x,y)
    use psb_base_mod, only : psb_dpk_, dzero, done
    real(psb_dpk_) ::  a2_box
    real(psb_dpk_), intent(in) :: x,y
    a2_box=done*epsilon
  end function a2_box
  function g_box(x,y)
    use psb_base_mod, only : psb_dpk_, dzero, done
    real(psb_dpk_) ::  g_box
    real(psb_dpk_), intent(in) :: x,y
    g_box = dzero
    if (x == done) then
      g_box = done
    else if (x == dzero) then
      g_box = done
    end if
  end function g_box
end module amg_d_pde2d_box_mod

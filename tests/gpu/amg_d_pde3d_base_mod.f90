module amg_d_pde3d_base_mod
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
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=done/sqrt(3.0_psb_dpk_)
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=done/sqrt(3.0_psb_dpk_)
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z
    b3=done/sqrt(3.0_psb_dpk_)
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z
    c=dzero
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a1
    real(psb_dpk_), intent(in) :: x,y,z
    a1=epsilon
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=epsilon
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=epsilon
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then
      g = done
    end if
  end function g
end module amg_d_pde3d_base_mod

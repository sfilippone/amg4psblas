module amg_s_pde3d_base_mod
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
  function b1(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y,z
    b1=sone/sqrt(3.0_psb_spk_)
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y,z
    b2=sone/sqrt(3.0_psb_spk_)
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) ::  b3
    real(psb_spk_), intent(in) :: x,y,z
    b3=sone/sqrt(3.0_psb_spk_)
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_spk_, sone
    real(psb_spk_) ::  c
    real(psb_spk_), intent(in) :: x,y,z
    c=szero
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a1
    real(psb_spk_), intent(in) :: x,y,z
    a1=epsilon
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a2
    real(psb_spk_), intent(in) :: x,y,z
    a2=epsilon
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_spk_
    real(psb_spk_) ::  a3
    real(psb_spk_), intent(in) :: x,y,z
    a3=epsilon
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    real(psb_spk_) ::  g
    real(psb_spk_), intent(in) :: x,y,z
    g = szero
    if (x == sone) then
      g = sone
    else if (x == szero) then
      g = sone
    end if
  end function g
end module amg_s_pde3d_base_mod

module amg_s_pde2d_base_mod
  use psb_base_mod, only : psb_spk_, szero, sone
  real(psb_spk_), save, private :: epsilon=sone/80
contains
  subroutine pde_set_parm(dat)
    real(psb_spk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm
  !
  ! functions parametrizing the differential equation
  !
  function b1(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y
    b1 = szero/1.414_psb_spk_
  end function b1
  function b2(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y
    b2 = szero/1.414_psb_spk_
  end function b2
  function c(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  c
    real(psb_spk_), intent(in) :: x,y
    c = szero
  end function c
  function a1(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  a1
    real(psb_spk_), intent(in) :: x,y
    a1=sone*epsilon
  end function a1
  function a2(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  a2
    real(psb_spk_), intent(in) :: x,y
    a2=sone*epsilon
  end function a2
  function g(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  g
    real(psb_spk_), intent(in) :: x,y
    g = szero
    if (x == sone) then
      g = sone
    else if (x == szero) then
      g = sone
    end if
  end function g
end module amg_s_pde2d_base_mod

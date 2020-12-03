module amg_s_pde2d_box_mod
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
  function b1_box(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) :: b1_box
    real(psb_spk_), intent(in) :: x,y
    b1_box = sone/1.414_psb_spk_
  end function b1_box
  function b2_box(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  b2_box
    real(psb_spk_), intent(in) :: x,y
    b2_box = sone/1.414_psb_spk_
  end function b2_box
  function c_box(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  c_box
    real(psb_spk_), intent(in) :: x,y
    c_box = szero
  end function c_box
  function a1_box(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  a1_box
    real(psb_spk_), intent(in) :: x,y
    a1_box=sone*epsilon
  end function a1_box
  function a2_box(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  a2_box
    real(psb_spk_), intent(in) :: x,y
    a2_box=sone*epsilon
  end function a2_box
  function g_box(x,y)
    use psb_base_mod, only : psb_spk_, szero, sone
    real(psb_spk_) ::  g_box
    real(psb_spk_), intent(in) :: x,y
    g_box = szero
    if (x == sone) then
      g_box = sone
    else if (x == szero) then
      g_box = sone
    end if
  end function g_box
end module amg_s_pde2d_box_mod

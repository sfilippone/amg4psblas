! C extensions for optional AMG parameters and diagnostics.
! The numerical operations are delegated to the existing native methods.
module amg_cbind_ext_mod
  use iso_c_binding
  use psb_base_cbind_mod
  use amg_prec_mod
  use amg_d_jac_smoother, only : amg_d_jac_smoother_type
  use amg_z_jac_smoother, only : amg_z_jac_smoother_type
  implicit none
  private

  integer, parameter :: amg_c_name_length = 80
  type, bind(c) :: amg_c_prec_info
    integer(psb_c_ipk_) :: nlevels, hierarchy_built, outer_sweeps
    real(c_double) :: operator_complexity, average_coarsening_ratio
    integer(c_int64_t) :: local_bytes
  end type amg_c_prec_info

  type, bind(c) :: amg_c_prec_level_info
    integer(psb_c_ipk_) :: cycle, coarse_matrix, sweeps_pre, sweeps_post
    integer(psb_c_ipk_) :: residual_stop_pre, residual_stop_post
    character(kind=c_char) :: smoother_pre(amg_c_name_length)
    character(kind=c_char) :: solver_pre(amg_c_name_length)
    character(kind=c_char) :: smoother_post(amg_c_name_length)
    character(kind=c_char) :: solver_post(amg_c_name_length)
  end type amg_c_prec_level_info

  public :: amg_c_prec_info, amg_c_prec_level_info
  public :: amg_c_dprecseti_opt, amg_c_dprecsetr_opt, amg_c_dprecsetc_opt
  public :: amg_c_dprec_get_info, amg_c_dprec_get_level_info
  public :: amg_c_zprecseti_opt, amg_c_zprecsetr_opt, amg_c_zprecsetc_opt
  public :: amg_c_zprec_get_info, amg_c_zprec_get_level_info

contains

  ! NULL leaves an allocatable optional argument unallocated, hence absent
  ! when forwarded to a non-allocatable OPTIONAL native dummy argument.
  subroutine read_string(cp, s, info)
    type(c_ptr), value :: cp
    character(len=:), allocatable, intent(out) :: s
    integer(psb_ipk_), intent(out) :: info
    character(kind=c_char), pointer :: p(:)
    integer :: i, n
    info = 0
    if (.not.c_associated(cp)) return
    call c_f_pointer(cp, p, [amg_c_name_length + 1])
    do n = 1, amg_c_name_length + 1
      if (p(n) == c_null_char) then
        allocate(character(len=n-1) :: s)
        do i = 1, n-1
          s(i:i) = p(i)
        end do
        return
      end if
    end do
    info = -1
  end subroutine read_string

  subroutine read_integer(cp, v)
    type(c_ptr), value :: cp
    integer(psb_ipk_), allocatable, intent(out) :: v
    integer(psb_c_ipk_), pointer :: p
    if (.not.c_associated(cp)) return
    call c_f_pointer(cp, p)
    allocate(v)
    v = int(p, kind=psb_ipk_)
  end subroutine read_integer

  subroutine read_options(cilev, cilmax, cpos, cidx, nlev, ilev, ilmax, pos, idx, info)
    type(c_ptr), value :: cilev, cilmax, cpos, cidx
    integer(psb_ipk_), intent(in) :: nlev
    integer(psb_ipk_), allocatable, intent(out) :: ilev, ilmax, idx
    character(len=:), allocatable, intent(out) :: pos
    integer(psb_ipk_), intent(out) :: info
    call read_integer(cilev, ilev)
    call read_integer(cilmax, ilmax)
    call read_integer(cidx, idx)
    call read_string(cpos, pos, info)
    if (info /= 0) return
    if (allocated(ilev)) then
      if ((ilev < 1).or.(ilev > nlev)) then
        info = -1
        return
      end if
      if (allocated(ilmax)) then
        if ((ilmax < ilev).or.(ilmax > nlev)) info = -1
      end if
    end if
  end subroutine read_options

  subroutine write_string(s, cp)
    character(len=*), intent(in) :: s
    character(kind=c_char), intent(out) :: cp(:)
    integer :: i
    cp = c_null_char
    do i = 1, min(len_trim(s), size(cp)-1)
      cp(i) = s(i:i)
    end do
  end subroutine write_string

  subroutine clear_info(out)
    type(amg_c_prec_info), intent(out) :: out
    out%nlevels = 0
    out%hierarchy_built = 0
    out%outer_sweeps = 0
    out%operator_complexity = 0
    out%average_coarsening_ratio = 0
    out%local_bytes = 0
  end subroutine clear_info

  subroutine clear_level_info(out)
    type(amg_c_prec_level_info), intent(out) :: out
    out%cycle = -1
    out%coarse_matrix = -1
    out%sweeps_pre = 0
    out%sweeps_post = 0
    out%residual_stop_pre = -1
    out%residual_stop_post = -1
    out%smoother_pre = c_null_char
    out%solver_pre = c_null_char
    out%smoother_post = c_null_char
    out%solver_post = c_null_char
  end subroutine clear_level_info

  function amg_c_dprecseti_opt(cph, what, val, cilev, cilmax, cpos, cidx) bind(c) result(res)
    type(c_ptr), value :: cph, what, cilev, cilmax, cpos, cidx
    integer(psb_c_ipk_), value :: val
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_dprec_type), pointer :: precp
    integer(psb_ipk_) :: info
    integer(psb_ipk_), allocatable :: ilev, ilmax, idx
    character(len=:), allocatable :: fwhat, pos

    res = -1
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if (.not.c_associated(what)) return
    call read_string(what, fwhat, info)
    if (info /= 0) return
    call read_options(cilev, cilmax, cpos, cidx, size(precp%precv, kind=psb_ipk_), &
         & ilev, ilmax, pos, idx, info)
    if (info /= 0) return
    call precp%set(fwhat, int(val, kind=psb_ipk_), info, ilev=ilev, ilmax=ilmax, pos=pos, idx=idx)
    res = int(info, kind=psb_c_ipk_)
  end function amg_c_dprecseti_opt

  function amg_c_dprecsetr_opt(cph, what, val, cilev, cilmax, cpos, cidx) bind(c) result(res)
    type(c_ptr), value :: cph, what, cilev, cilmax, cpos, cidx
    real(c_double), value :: val
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_dprec_type), pointer :: precp
    integer(psb_ipk_) :: info
    integer(psb_ipk_), allocatable :: ilev, ilmax, idx
    character(len=:), allocatable :: fwhat, pos

    res = -1
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if (.not.c_associated(what)) return
    call read_string(what, fwhat, info)
    if (info /= 0) return
    call read_options(cilev, cilmax, cpos, cidx, size(precp%precv, kind=psb_ipk_), &
         & ilev, ilmax, pos, idx, info)
    if (info /= 0) return
    call precp%set(fwhat, real(val, kind=psb_dpk_), info, ilev=ilev, ilmax=ilmax, pos=pos, idx=idx)
    res = int(info, kind=psb_c_ipk_)
  end function amg_c_dprecsetr_opt

  function amg_c_dprecsetc_opt(cph, what, val, cilev, cilmax, cpos, cidx) bind(c) result(res)
    type(c_ptr), value :: cph, what, cilev, cilmax, cpos, cidx
    type(c_ptr), value :: val
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_dprec_type), pointer :: precp
    integer(psb_ipk_) :: info
    integer(psb_ipk_), allocatable :: ilev, ilmax, idx
    character(len=:), allocatable :: fwhat, pos, fval

    res = -1
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if (.not.c_associated(what)) return
    call read_string(what, fwhat, info)
    if (info /= 0) return
    if (.not.c_associated(val)) return
    call read_string(val, fval, info)
    if (info /= 0) return
    call read_options(cilev, cilmax, cpos, cidx, size(precp%precv, kind=psb_ipk_), &
         & ilev, ilmax, pos, idx, info)
    if (info /= 0) return
    call precp%set(fwhat, fval, info, ilev=ilev, ilmax=ilmax, pos=pos, idx=idx)
    res = int(info, kind=psb_c_ipk_)
  end function amg_c_dprecsetc_opt

  function amg_c_dprec_get_info(cph, cout) bind(c) result(res)
    type(c_ptr), value :: cph, cout
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_dprec_type), pointer :: precp
    type(amg_c_prec_info), pointer :: out

    res = -1
    if (.not.c_associated(cout)) return
    call c_f_pointer(cout, out)
    call clear_info(out)
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    out%outer_sweeps = int(precp%outer_sweeps, kind=psb_c_ipk_)
    out%local_bytes = int(precp%sizeof(global=.false.), kind=c_int64_t)
    if (allocated(precp%precv)) then
      out%nlevels = size(precp%precv, kind=psb_c_ipk_)
      if (out%nlevels > 0) then
        if (associated(precp%precv(1)%base_a)) then
          out%hierarchy_built = 1
          out%nlevels = int(precp%get_nlevs(), kind=psb_c_ipk_)
          out%operator_complexity = real(precp%get_complexity(), kind=c_double)
          out%average_coarsening_ratio = real(precp%get_avg_cr(), kind=c_double)
        end if
      end if
    end if
    res = 0
  end function amg_c_dprec_get_info

  function amg_c_dprec_get_level_info(cph, ilev, cout) bind(c) result(res)
    type(c_ptr), value :: cph, cout
    integer(psb_c_ipk_), value :: ilev
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_dprec_type), pointer :: precp
    type(amg_c_prec_level_info), pointer :: out

    res = -1
    if (.not.c_associated(cout)) return
    call c_f_pointer(cout, out)
    call clear_level_info(out)
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if ((ilev < 1).or.(ilev > size(precp%precv))) return
    if (associated(precp%precv(1)%base_a)) then
      if (ilev > precp%get_nlevs()) return
    end if
    associate(lv => precp%precv(ilev))
      out%cycle = int(lv%parms%ml_cycle, kind=psb_c_ipk_)
      out%coarse_matrix = int(lv%parms%coarse_mat, kind=psb_c_ipk_)
      out%sweeps_pre = int(lv%parms%sweeps_pre, kind=psb_c_ipk_)
      out%sweeps_post = int(lv%parms%sweeps_post, kind=psb_c_ipk_)
      if (allocated(lv%sm)) call d_smoother_info(lv%sm, &
           & out%smoother_pre, out%solver_pre, out%residual_stop_pre)
      if (allocated(lv%sm2a)) then
        call d_smoother_info(lv%sm2a, &
             & out%smoother_post, out%solver_post, out%residual_stop_post)
      else if (allocated(lv%sm)) then
        ! With no separate post-smoother, native setup aliases it to sm.
        call d_smoother_info(lv%sm, &
             & out%smoother_post, out%solver_post, out%residual_stop_post)
      end if
    end associate
    res = 0
  end function amg_c_dprec_get_level_info

  subroutine d_smoother_info(sm, name, solver, stop)
    class(amg_d_base_smoother_type), intent(in) :: sm
    character(kind=c_char), intent(out) :: name(:), solver(:)
    integer(psb_c_ipk_), intent(out) :: stop
    call write_string(sm%get_fmt(), name)
    solver = c_null_char
    if (allocated(sm%sv)) call write_string(sm%sv%get_fmt(), solver)
    stop = -1
    select type(sm)
    class is(amg_d_jac_smoother_type)
      stop = 0
      if (sm%checkres) stop = 1
    end select
  end subroutine d_smoother_info

  function amg_c_zprecseti_opt(cph, what, val, cilev, cilmax, cpos, cidx) bind(c) result(res)
    type(c_ptr), value :: cph, what, cilev, cilmax, cpos, cidx
    integer(psb_c_ipk_), value :: val
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_zprec_type), pointer :: precp
    integer(psb_ipk_) :: info
    integer(psb_ipk_), allocatable :: ilev, ilmax, idx
    character(len=:), allocatable :: fwhat, pos

    res = -1
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if (.not.c_associated(what)) return
    call read_string(what, fwhat, info)
    if (info /= 0) return
    call read_options(cilev, cilmax, cpos, cidx, size(precp%precv, kind=psb_ipk_), &
         & ilev, ilmax, pos, idx, info)
    if (info /= 0) return
    call precp%set(fwhat, int(val, kind=psb_ipk_), info, ilev=ilev, ilmax=ilmax, pos=pos, idx=idx)
    res = int(info, kind=psb_c_ipk_)
  end function amg_c_zprecseti_opt

  function amg_c_zprecsetr_opt(cph, what, val, cilev, cilmax, cpos, cidx) bind(c) result(res)
    type(c_ptr), value :: cph, what, cilev, cilmax, cpos, cidx
    real(c_double), value :: val
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_zprec_type), pointer :: precp
    integer(psb_ipk_) :: info
    integer(psb_ipk_), allocatable :: ilev, ilmax, idx
    character(len=:), allocatable :: fwhat, pos

    res = -1
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if (.not.c_associated(what)) return
    call read_string(what, fwhat, info)
    if (info /= 0) return
    call read_options(cilev, cilmax, cpos, cidx, size(precp%precv, kind=psb_ipk_), &
         & ilev, ilmax, pos, idx, info)
    if (info /= 0) return
    call precp%set(fwhat, real(val, kind=psb_dpk_), info, ilev=ilev, ilmax=ilmax, pos=pos, idx=idx)
    res = int(info, kind=psb_c_ipk_)
  end function amg_c_zprecsetr_opt

  function amg_c_zprecsetc_opt(cph, what, val, cilev, cilmax, cpos, cidx) bind(c) result(res)
    type(c_ptr), value :: cph, what, cilev, cilmax, cpos, cidx
    type(c_ptr), value :: val
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_zprec_type), pointer :: precp
    integer(psb_ipk_) :: info
    integer(psb_ipk_), allocatable :: ilev, ilmax, idx
    character(len=:), allocatable :: fwhat, pos, fval

    res = -1
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if (.not.c_associated(what)) return
    call read_string(what, fwhat, info)
    if (info /= 0) return
    if (.not.c_associated(val)) return
    call read_string(val, fval, info)
    if (info /= 0) return
    call read_options(cilev, cilmax, cpos, cidx, size(precp%precv, kind=psb_ipk_), &
         & ilev, ilmax, pos, idx, info)
    if (info /= 0) return
    call precp%set(fwhat, fval, info, ilev=ilev, ilmax=ilmax, pos=pos, idx=idx)
    res = int(info, kind=psb_c_ipk_)
  end function amg_c_zprecsetc_opt

  function amg_c_zprec_get_info(cph, cout) bind(c) result(res)
    type(c_ptr), value :: cph, cout
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_zprec_type), pointer :: precp
    type(amg_c_prec_info), pointer :: out

    res = -1
    if (.not.c_associated(cout)) return
    call c_f_pointer(cout, out)
    call clear_info(out)
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    out%outer_sweeps = int(precp%outer_sweeps, kind=psb_c_ipk_)
    out%local_bytes = int(precp%sizeof(global=.false.), kind=c_int64_t)
    if (allocated(precp%precv)) then
      out%nlevels = size(precp%precv, kind=psb_c_ipk_)
      if (out%nlevels > 0) then
        if (associated(precp%precv(1)%base_a)) then
          out%hierarchy_built = 1
          out%nlevels = int(precp%get_nlevs(), kind=psb_c_ipk_)
          out%operator_complexity = real(precp%get_complexity(), kind=c_double)
          out%average_coarsening_ratio = real(precp%get_avg_cr(), kind=c_double)
        end if
      end if
    end if
    res = 0
  end function amg_c_zprec_get_info

  function amg_c_zprec_get_level_info(cph, ilev, cout) bind(c) result(res)
    type(c_ptr), value :: cph, cout
    integer(psb_c_ipk_), value :: ilev
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), pointer :: ph
    type(amg_zprec_type), pointer :: precp
    type(amg_c_prec_level_info), pointer :: out

    res = -1
    if (.not.c_associated(cout)) return
    call c_f_pointer(cout, out)
    call clear_level_info(out)
    if (.not.c_associated(cph)) return
    call c_f_pointer(cph, ph)
    if (.not.c_associated(ph%item)) return
    call c_f_pointer(ph%item, precp)
    if (.not.allocated(precp%precv)) return
    if ((ilev < 1).or.(ilev > size(precp%precv))) return
    if (associated(precp%precv(1)%base_a)) then
      if (ilev > precp%get_nlevs()) return
    end if
    associate(lv => precp%precv(ilev))
      out%cycle = int(lv%parms%ml_cycle, kind=psb_c_ipk_)
      out%coarse_matrix = int(lv%parms%coarse_mat, kind=psb_c_ipk_)
      out%sweeps_pre = int(lv%parms%sweeps_pre, kind=psb_c_ipk_)
      out%sweeps_post = int(lv%parms%sweeps_post, kind=psb_c_ipk_)
      if (allocated(lv%sm)) call z_smoother_info(lv%sm, &
           & out%smoother_pre, out%solver_pre, out%residual_stop_pre)
      if (allocated(lv%sm2a)) then
        call z_smoother_info(lv%sm2a, &
             & out%smoother_post, out%solver_post, out%residual_stop_post)
      else if (allocated(lv%sm)) then
        ! With no separate post-smoother, native setup aliases it to sm.
        call z_smoother_info(lv%sm, &
             & out%smoother_post, out%solver_post, out%residual_stop_post)
      end if
    end associate
    res = 0
  end function amg_c_zprec_get_level_info

  subroutine z_smoother_info(sm, name, solver, stop)
    class(amg_z_base_smoother_type), intent(in) :: sm
    character(kind=c_char), intent(out) :: name(:), solver(:)
    integer(psb_c_ipk_), intent(out) :: stop
    call write_string(sm%get_fmt(), name)
    solver = c_null_char
    if (allocated(sm%sv)) call write_string(sm%sv%get_fmt(), solver)
    stop = -1
    select type(sm)
    class is(amg_z_jac_smoother_type)
      stop = 0
      if (sm%checkres) stop = 1
    end select
  end subroutine z_smoother_info

end module amg_cbind_ext_mod

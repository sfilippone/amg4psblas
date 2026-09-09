module psb_mixed_support_mod
  use psb_base_mod

contains

  subroutine psb_d2s_cscnv(a,b,info,type,mold)
    implicit none 
    class(psb_dspmat_type), intent(in)    :: a
    class(psb_sspmat_type), intent(out)   :: b
    integer(psb_ipk_), intent(out)                   :: info
    character(len=*), optional, intent(in) :: type
    class(psb_s_base_sparse_mat), intent(in), optional :: mold

    type(psb_d_coo_sparse_mat) :: dcoo
    type(psb_s_coo_sparse_mat) :: scoo
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_erractionsave(err_act)

    call a%cp_to(dcoo)
    call coo_d2s(dcoo,scoo,info)

    if (present(mold)) then

      allocate(b%a, mold=mold,stat=info)

    else if (present(type)) then

      select case (psb_toupper(type))
      case ('CSR')
        allocate(psb_s_csr_sparse_mat :: b%a, stat=info)
      case ('COO')
        allocate(psb_s_coo_sparse_mat :: b%a, stat=info)
      case ('CSC')
        allocate(psb_s_csc_sparse_mat :: b%a, stat=info)
      case default
        info = psb_err_format_unknown_
        call psb_errpush(info,name,a_err=type)
        goto 9999
      end select
    end if
    call b%a%mv_from_coo(scoo,info)
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  contains
    subroutine coo_d2s(dcoo,scoo,info)
      type(psb_d_coo_sparse_mat) :: dcoo
      type(psb_s_coo_sparse_mat) :: scoo
      integer(psb_ipk_), intent(out)                   :: info

      integer(psb_ipk_) :: i,j,nr,nc,nz
      info = psb_success_
      nr = dcoo%get_nrows()
      nc = dcoo%get_ncols()
      nz = dcoo%get_nzeros()
      call scoo%allocate(nr,nc,nz)
      do i=1,nz
        scoo%ia(i) = dcoo%ia(i)
        scoo%ja(i) = dcoo%ja(i)
        scoo%val(i) = dcoo%val(i)
      end do
      call scoo%set_nzeros(nz)
      return
    end subroutine coo_d2s

  end subroutine psb_d2s_cscnv

  subroutine psb_s2d_cscnv(a,b,info,type,mold)
    implicit none 
    class(psb_sspmat_type), intent(in)    :: a
    class(psb_dspmat_type), intent(out)   :: b
    integer(psb_ipk_), intent(out)                   :: info
    character(len=*), optional, intent(in) :: type
    class(psb_d_base_sparse_mat), intent(in), optional :: mold

    type(psb_d_coo_sparse_mat) :: dcoo
    type(psb_s_coo_sparse_mat) :: scoo
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_erractionsave(err_act)

    call a%cp_to(scoo)
    call coo_s2d(scoo,dcoo,info)

    if (present(mold)) then

      allocate(b%a, mold=mold,stat=info)

    else if (present(type)) then

      select case (psb_toupper(type))
      case ('CSR')
        allocate(psb_d_csr_sparse_mat :: b%a, stat=info)
      case ('COO')
        allocate(psb_d_coo_sparse_mat :: b%a, stat=info)
      case ('CSC')
        allocate(psb_d_csc_sparse_mat :: b%a, stat=info)
      case default
        info = psb_err_format_unknown_
        call psb_errpush(info,name,a_err=type)
        goto 9999
      end select
    end if
    call b%a%mv_from_coo(dcoo,info)
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  contains
    subroutine coo_s2d(scoo,dcoo,info)
      type(psb_d_coo_sparse_mat) :: dcoo
      type(psb_s_coo_sparse_mat) :: scoo
      integer(psb_ipk_), intent(out)                   :: info

      integer(psb_ipk_) :: i,j,nr,nc,nz
      info = psb_success_
      nr = scoo%get_nrows()
      nc = scoo%get_ncols()
      nz = scoo%get_nzeros()
      call scoo%allocate(nr,nc,nz)
      do i=1,nz
        dcoo%ia(i) = scoo%ia(i)
        dcoo%ja(i) = scoo%ja(i)
        dcoo%val(i) = scoo%val(i)
      end do
      call dcoo%set_nzeros(nz)
      return
    end subroutine coo_s2d

  end subroutine psb_s2d_cscnv

  subroutine psb_d2s_vect(dv,sv,info,mold)
    class(psb_d_vect_type), intent(inout) :: dv
    class(psb_s_vect_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)        :: info    
    class(psb_s_base_vect_type), intent(in), optional :: mold

    real(psb_spk_), allocatable :: xv(:)

    info = psb_success_
    xv = dv%get_vect()
    call sv%bld(xv,mold=mold)
  end subroutine psb_d2s_vect

  subroutine psb_s2d_vect(sv,dv,info,mold)
    class(psb_d_vect_type), intent(inout) :: dv
    class(psb_s_vect_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)        :: info    
    class(psb_d_base_vect_type), intent(in), optional :: mold

    real(psb_dpk_), allocatable :: xv(:)

    info = psb_success_
    xv = sv%get_vect()
    call dv%bld(xv,mold=mold)
  end subroutine psb_s2d_vect

end module psb_mixed_support_mod

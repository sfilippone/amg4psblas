!
!
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.7)
!
!
subroutine amg_d_richards_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,&
     & sweeps,work,wv,info,init,initu)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_apply_vect
  implicit none

  type(psb_desc_type), intent(in)                    :: desc_data
  class(amg_d_richards_smoother_type), intent(inout) :: sm
  type(psb_d_vect_type), intent(inout)               :: x
  type(psb_d_vect_type), intent(inout)               :: y
  real(psb_dpk_), intent(in)                         :: alpha, beta
  character(len=1), intent(in)                       :: trans
  integer(psb_ipk_), intent(in)                      :: sweeps
  real(psb_dpk_), target, intent(inout)              :: work(:)
  type(psb_d_vect_type), intent(inout)               :: wv(:)
  integer(psb_ipk_), intent(out)                     :: info
  character, intent(in), optional                    :: init
  type(psb_d_vect_type), intent(inout), optional     :: initu

  integer(psb_ipk_)   :: n_col, err_act, i
  character           :: trans_, init_
  real(psb_dpk_), pointer :: aux(:)
  character(len=32)   :: name='d_richards_smoother_apply_v'

  call psb_erractionsave(err_act)
  info = psb_success_

  trans_ = psb_toupper(trans)
  if ((trans_ /= 'N').and.(trans_ /= 'T').and.(trans_ /= 'C')) then
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end if

  if (sweeps < 0) then
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/itwo,sweeps,izero,izero,izero/))
    goto 9999
  end if

  if (.not.allocated(sm%sv)) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (.not.associated(sm%pa)) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name,a_err='matrix pointer not associated')
    goto 9999
  end if

  if (size(wv) < 3) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='workspace vectors too small')
    goto 9999
  end if

  init_ = 'Z'
  if (present(init)) init_ = psb_toupper(init)

  n_col = desc_data%get_local_cols()
  if (4*n_col <= size(work)) then
    aux => work(:)
  else
    allocate(aux(4*n_col),stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/4*n_col,izero,izero,izero,izero/),&
           & a_err='real(psb_dpk_)')
      goto 9999
    end if
  end if

  associate(tx => wv(1), ty => wv(2), tz => wv(3))

    select case(init_)
    case('Z')
      call psb_geaxpby(dzero,x,dzero,ty,desc_data,info)
    case('Y')
      call psb_geaxpby(done,y,dzero,ty,desc_data,info)
    case('U')
      if (.not.present(initu)) then
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='missing initu to smoother_apply_vect')
        goto 9999
      end if
      call psb_geaxpby(done,initu,dzero,ty,desc_data,info)
    case default
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='wrong init to smoother_apply_vect')
      goto 9999
    end select

    do i = 1, sweeps
      call psb_geaxpby(done,x,dzero,tx,desc_data,info)
      call psb_spmm(-done,sm%pa,ty,done,tx,desc_data,info,work=aux,trans=trans_)
      if (info /= psb_success_) exit

      call sm%sv%apply(done,tx,dzero,tz,desc_data,trans_,aux,wv(4:),info,init='Y')
      if (info /= psb_success_) exit

      call psb_geaxpby(sm%omega,tz,done,ty,desc_data,info)
      if (info /= psb_success_) exit
    end do

    if (info == psb_success_) call psb_geaxpby(alpha,ty,beta,y,desc_data,info)

  end associate

  if (.not.(4*n_col <= size(work))) deallocate(aux)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_apply_vect


subroutine amg_d_richards_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,&
     & sweeps,work,info,init,initu)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_apply
  implicit none

  type(psb_desc_type), intent(in)                    :: desc_data
  class(amg_d_richards_smoother_type), intent(inout) :: sm
  real(psb_dpk_), intent(inout)                      :: x(:)
  real(psb_dpk_), intent(inout)                      :: y(:)
  real(psb_dpk_), intent(in)                         :: alpha, beta
  character(len=1), intent(in)                       :: trans
  integer(psb_ipk_), intent(in)                      :: sweeps
  real(psb_dpk_), target, intent(inout)              :: work(:)
  integer(psb_ipk_), intent(out)                     :: info
  character, intent(in), optional                    :: init
  real(psb_dpk_), intent(inout), optional            :: initu(:)

  integer(psb_ipk_)      :: n_col, err_act, i
  character              :: trans_, init_
  real(psb_dpk_), pointer :: aux(:)
  real(psb_dpk_), allocatable :: tx(:), ty(:), tz(:)
  character(len=30)      :: name='d_richards_smoother_apply'

  call psb_erractionsave(err_act)
  info = psb_success_

  trans_ = psb_toupper(trans)
  if ((trans_ /= 'N').and.(trans_ /= 'T').and.(trans_ /= 'C')) then
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end if

  if (sweeps < 0) then
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/itwo,sweeps,izero,izero,izero/))
    goto 9999
  end if

  if (.not.allocated(sm%sv)) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (.not.associated(sm%pa)) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name,a_err='matrix pointer not associated')
    goto 9999
  end if

  n_col = desc_data%get_local_cols()
  if (4*n_col <= size(work)) then
    aux => work(:)
  else
    allocate(aux(4*n_col),stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/4*n_col,izero,izero,izero,izero/),&
           & a_err='real(psb_dpk_)')
      goto 9999
    end if
  end if

  allocate(tx(size(x)),ty(size(y)),tz(size(y)),stat=info)
  if (info /= 0) then
    info = psb_err_alloc_request_
    call psb_errpush(info,name,a_err='alloc tx/ty/tz')
    goto 9999
  end if

  init_ = 'Z'
  if (present(init)) init_ = psb_toupper(init)

  select case(init_)
  case('Z')
    ty(:) = dzero
  case('Y')
    ty(:) = y(:)
  case('U')
    if (.not.present(initu)) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='missing initu to smoother_apply')
      goto 9999
    end if
    ty(:) = initu(:)
  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='wrong init to smoother_apply')
    goto 9999
  end select

  do i = 1, sweeps
    tx(:) = x(:)
    call psb_spmm(-done,sm%pa,ty,done,tx,desc_data,info,work=aux,trans=trans_)
    if (info /= psb_success_) exit

    call sm%sv%apply(done,tx,dzero,tz,desc_data,trans_,aux,info,init='Y')
    if (info /= psb_success_) exit

    ty(:) = ty(:) + sm%omega*tz(:)
  end do

  if (info == psb_success_) y(:) = beta*y(:) + alpha*ty(:)

  deallocate(tx,ty,tz,stat=info)
  if (.not.(4*n_col <= size(work))) deallocate(aux)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_apply


subroutine amg_d_richards_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_bld
  implicit none

  type(psb_dspmat_type), intent(inout), target         :: a
  type(psb_desc_type), intent(inout)                   :: desc_a
  class(amg_d_richards_smoother_type), intent(inout)   :: sm
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_d_base_sparse_mat), intent(in), optional   :: amold
  class(psb_d_base_vect_type), intent(in), optional    :: vmold
  class(psb_i_base_vect_type), intent(in), optional    :: imold

  integer(psb_ipk_) :: err_act
  character(len=28) :: name='d_richards_smoother_bld'

  info = psb_success_
  call psb_erractionsave(err_act)

  sm%pa => a
  if (sm%omega == dzero) sm%omega = done

  if (.not.allocated(sm%sv)) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call sm%sv%build(a,desc_a,info,amold=amold,vmold=vmold,imold=imold)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='solver build')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_bld


subroutine amg_d_richards_smoother_cnv(sm,info,amold,vmold,imold)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_cnv
  implicit none

  class(amg_d_richards_smoother_type), intent(inout)   :: sm
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_d_base_sparse_mat), intent(in), optional   :: amold
  class(psb_d_base_vect_type), intent(in), optional    :: vmold
  class(psb_i_base_vect_type), intent(in), optional    :: imold

  integer(psb_ipk_) :: err_act
  character(len=28) :: name='d_richards_smoother_cnv'

  info = psb_success_
  call psb_erractionsave(err_act)

  if (allocated(sm%sv)) call sm%sv%cnv(info,amold=amold,vmold=vmold,imold=imold)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='solver cnv')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_cnv


subroutine amg_d_richards_smoother_dmp(sm,desc,level,info,prefix,head,smoother,solver,global_num)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_dmp
  implicit none

  class(amg_d_richards_smoother_type), intent(in)      :: sm
  type(psb_desc_type), intent(in)                      :: desc
  integer(psb_ipk_), intent(in)                        :: level
  integer(psb_ipk_), intent(out)                       :: info
  character(len=*), intent(in), optional               :: prefix, head
  logical, intent(in), optional                        :: smoother, solver, global_num

  info = psb_success_
  if (allocated(sm%sv)) call sm%sv%dump(desc,level,info,solver=solver,prefix=prefix,global_num=global_num)

end subroutine amg_d_richards_smoother_dmp


subroutine amg_d_richards_smoother_clone(sm,smout,info)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_clone
  implicit none

  class(amg_d_richards_smoother_type), intent(inout)               :: sm
  class(amg_d_base_smoother_type), allocatable, intent(inout)      :: smout
  integer(psb_ipk_), intent(out)                                   :: info

  integer(psb_ipk_) :: err_act

  info = psb_success_
  call psb_erractionsave(err_act)

  if (allocated(smout)) then
    call smout%free(info)
    if (info == psb_success_) deallocate(smout, stat=info)
  end if
  if (info == psb_success_) allocate(amg_d_richards_smoother_type :: smout, stat=info)
  if (info /= 0) then
    info = psb_err_alloc_dealloc_
    goto 9999
  end if

  select type(smo => smout)
  type is (amg_d_richards_smoother_type)
    smo%global_nnz_tot = sm%global_nnz_tot
    smo%checkres = sm%checkres
    smo%printres = sm%printres
    smo%checkiter = sm%checkiter
    smo%printiter = sm%printiter
    smo%tol = sm%tol
    smo%omega = sm%omega
    smo%pa => sm%pa
    if (allocated(sm%sv)) then
      allocate(smo%sv,mold=sm%sv,stat=info)
      if (info == psb_success_) call sm%sv%clone(smo%sv,info)
    end if
  class default
    info = psb_err_internal_error_
  end select

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_clone


subroutine amg_d_richards_smoother_clone_settings(sm,smout,info)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_clone_settings
  implicit none

  class(amg_d_richards_smoother_type), intent(inout) :: sm
  class(amg_d_base_smoother_type), intent(inout)     :: smout
  integer(psb_ipk_), intent(out)                     :: info

  integer(psb_ipk_) :: err_act
  character(len=36) :: name='d_richards_smoother_clone_settings'

  call psb_erractionsave(err_act)
  info = psb_success_

  select type(smout)
  class is(amg_d_richards_smoother_type)
    smout%pa => null()
    smout%global_nnz_tot = 0
    smout%checkres = sm%checkres
    smout%printres = sm%printres
    smout%checkiter = sm%checkiter
    smout%printiter = sm%printiter
    smout%tol = sm%tol
    smout%omega = sm%omega

    if (allocated(smout%sv)) then
      if (.not.same_type_as(sm%sv,smout%sv)) then
        call smout%sv%free(info)
        if (info == 0) deallocate(smout%sv,stat=info)
      end if
    end if

    if (info == 0) then
      if (allocated(smout%sv)) then
        if (same_type_as(sm%sv,smout%sv)) then
          call sm%sv%clone_settings(smout%sv,info)
        else
          info = psb_err_internal_error_
        end if
      else
        allocate(smout%sv,mold=sm%sv,stat=info)
        if (info == 0) call sm%sv%clone_settings(smout%sv,info)
      end if
    end if
  class default
    info = psb_err_internal_error_
  end select

  if (info /= 0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_clone_settings


subroutine amg_d_richards_smoother_clear_data(sm,info)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_clear_data
  implicit none

  class(amg_d_richards_smoother_type), intent(inout) :: sm
  integer(psb_ipk_), intent(out)                     :: info

  integer(psb_ipk_) :: err_act
  character(len=31) :: name='d_richards_smoother_clear_data'

  call psb_erractionsave(err_act)

  info = psb_success_
  sm%pa => null()
  sm%global_nnz_tot = 0
  if ((info == 0).and.allocated(sm%sv)) call sm%sv%clear_data(info)
  if (info /= 0) then
    info = psb_err_internal_error_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_clear_data


subroutine amg_d_richards_smoother_descr(sm,info,iout,coarse,prefix)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_descr
  implicit none

  class(amg_d_richards_smoother_type), intent(in) :: sm
  integer(psb_ipk_), intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional         :: iout
  logical, intent(in), optional                   :: coarse
  character(len=*), intent(in), optional          :: prefix

  integer(psb_ipk_)   :: iout_
  character(1024)     :: prefix_

  info = psb_success_
  iout_ = psb_out_unit
  if (present(iout)) iout_ = iout
  prefix_ = ''
  if (present(prefix)) prefix_ = prefix

  write(iout_,*) trim(prefix_), '  Richardson smoother'
  write(iout_,*) trim(prefix_), '       Relaxation omega: ', sm%omega
  if (allocated(sm%sv)) then
    write(iout_,*) trim(prefix_), '       Preconditioner details:'
    call sm%sv%descr(info,iout_,coarse=coarse,prefix=prefix)
  end if

end subroutine amg_d_richards_smoother_descr


subroutine amg_d_richards_smoother_cseti(sm,what,val,info,idx)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_cseti
  implicit none

  class(amg_d_richards_smoother_type), intent(inout) :: sm
  character(len=*), intent(in)                       :: what
  integer(psb_ipk_), intent(in)                      :: val
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_), intent(in), optional            :: idx

  integer(psb_ipk_) :: err_act

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(what))
  case('SMOOTHER_RESIDUAL')
    sm%checkiter = val
  case('SMOOTHER_ITRACE')
    sm%printiter = val
  case default
    call sm%amg_d_base_smoother_type%set(what,val,info,idx=idx)
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_cseti


subroutine amg_d_richards_smoother_csetc(sm,what,val,info,idx)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_csetc
  implicit none

  class(amg_d_richards_smoother_type), intent(inout) :: sm
  character(len=*), intent(in)                       :: what
  character(len=*), intent(in)                       :: val
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_), intent(in), optional            :: idx

  integer(psb_ipk_) :: err_act
  character(len=30) :: name='d_richards_smoother_csetc'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(trim(what)))
  case('SMOOTHER_STOP')
    select case(psb_toupper(trim(val)))
    case('T','TRUE')
      sm%checkres = .true.
    case('F','FALSE')
      sm%checkres = .false.
    end select
  case('SMOOTHER_TRACE')
    select case(psb_toupper(trim(val)))
    case('T','TRUE')
      sm%printres = .true.
    case('F','FALSE')
      sm%printres = .false.
    end select
  case default
    call sm%amg_d_base_smoother_type%set(what,val,info,idx=idx)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_csetc


subroutine amg_d_richards_smoother_csetr(sm,what,val,info,idx)

  use psb_base_mod
  use amg_d_richards_smoother, amg_protect_name => amg_d_richards_smoother_csetr
  implicit none

  class(amg_d_richards_smoother_type), intent(inout) :: sm
  character(len=*), intent(in)                       :: what
  real(psb_dpk_), intent(in)                         :: val
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_), intent(in), optional            :: idx

  integer(psb_ipk_) :: err_act

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(what))
  case('SMOOTHER_STOPTOL')
    sm%tol = val
  case('RICHARDSON_OMEGA','SMOOTHER_OMEGA','OMEGA')
    sm%omega = val
  case default
    call sm%amg_d_base_smoother_type%set(what,val,info,idx=idx)
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_d_richards_smoother_csetr

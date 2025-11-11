#if !defined(PSB_IPK8)
subroutine amg_z_slu_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

  use psb_base_mod
  use amg_z_slu_solver, amg_protect_name => amg_z_slu_solver_bld
  Implicit None

  ! Arguments
  type(psb_zspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(inout)                  :: desc_a 
  class(amg_z_slu_solver_type), intent(inout)         :: sv
  integer(psb_ipk_), intent(out)                      :: info
  type(psb_zspmat_type), intent(in), target, optional :: b
  class(psb_z_base_sparse_mat), intent(in), optional  :: amold
  class(psb_z_base_vect_type), intent(in), optional   :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  ! Local variables
  type(psb_zspmat_type) :: atmp
  type(psb_z_csc_sparse_mat) :: acsc
  type(psb_z_coo_sparse_mat) :: acoo
  integer :: n_row,n_col, nrow_a, nztota
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)  :: np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='s_slu_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()
  n_col  = desc_a%get_local_cols()


  call a%cscnv(atmp,info,type='coo')
  call psb_rwextd(n_row,atmp,info,b=b) 
  call atmp%cscnv(info,type='coo',dupl=psb_dupl_add_)
  nrow_a = atmp%get_nrows()
  call atmp%a%csclip(acoo,info,jmax=nrow_a)
  call acsc%mv_from_coo(acoo,info)
  nztota = acsc%get_nzeros()
  ! Fix the entries to call C-base SuperLU
  acsc%ia(:)  = acsc%ia(:)  - 1
  acsc%icp(:) = acsc%icp(:) - 1
  info = amg_zslu_fact(nrow_a,nztota,acsc%val,&
       & acsc%icp,acsc%ia,sv%lufactors)

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='amg_zslu_fact'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call acsc%free()
  call atmp%free()

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine amg_z_slu_solver_bld

subroutine amg_z_slu_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
     & trans,work,wv,info,init,initu)
  use psb_base_mod
  use amg_z_slu_solver, amg_protect_name => amg_z_slu_solver_apply_vect

  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(amg_z_slu_solver_type), intent(inout) :: sv
  type(psb_z_vect_type),intent(inout)  :: x
  type(psb_z_vect_type),intent(inout)  :: y
  complex(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)           :: trans
  complex(psb_dpk_),target, intent(inout) :: work(:)
  type(psb_z_vect_type),intent(inout) :: wv(:)
  integer(psb_ipk_), intent(out)                :: info
  character, intent(in), optional                :: init
  type(psb_z_vect_type),intent(inout), optional   :: initu

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='s_slu_solver_apply_vect'

  call psb_erractionsave(err_act)

  info = psb_success_
  !
  ! For non-iterative solvers, init and initu are ignored.
  !

  call x%v%sync()
  call y%v%sync()
  call sv%apply(alpha,x%v%v,beta,y%v%v,desc_data,trans,work,info)
  call y%v%set_host()
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_z_slu_solver_apply_vect

subroutine amg_z_slu_solver_apply(alpha,sv,x,beta,y,desc_data,&
     & trans,work,info,init,initu)
  use psb_base_mod
  use amg_z_slu_solver, amg_protect_name => amg_z_slu_solver_apply
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(amg_z_slu_solver_type), intent(inout) :: sv
  complex(psb_dpk_),intent(inout)         :: x(:)
  complex(psb_dpk_),intent(inout)         :: y(:)
  complex(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)          :: trans
  complex(psb_dpk_),target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(out)        :: info
  character, intent(in), optional       :: init
  complex(psb_dpk_),intent(inout), optional :: initu(:)

  integer(psb_ipk_)  :: n_row,n_col
  complex(psb_dpk_), pointer :: ww(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='s_slu_solver_apply'

  call psb_erractionsave(err_act)

  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select
  !
  ! For non-iterative solvers, init and initu are ignored.
  !

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (n_col <= size(work)) then 
    ww => work(1:n_col)
  else
    allocate(ww(n_col),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/n_col,0,0,0,0/),&
           & a_err='complex(psb_dpk_)')
      goto 9999      
    end if
  endif

  ww(1:n_row) = x(1:n_row)
  select case(trans_)
  case('N')
    info = amg_zslu_solve(0,n_row,1,ww,n_row,sv%lufactors)
  case('T')
    info = amg_zslu_solve(1,n_row,1,ww,n_row,sv%lufactors)
  case('C')
    info = amg_zslu_solve(2,n_row,1,ww,n_row,sv%lufactors)
  case default
    call psb_errpush(psb_err_internal_error_, &
         & name,a_err='Invalid TRANS in ILU subsolve')
    goto 9999
  end select

  if (info == psb_success_) &
       & call psb_geaxpby(alpha,ww,beta,y,desc_data,info)


  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,& 
         & name,a_err='Error in subsolve')
    goto 9999
  endif

  if (n_col > size(work)) then 
    deallocate(ww)
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_z_slu_solver_apply
#endif

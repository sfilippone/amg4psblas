!  
!   
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2020 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Fabio Durastante        
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
!
! File: amg_c_slu_solver_mod.f90
!
! Module: amg_c_slu_solver_mod
!
!  This module defines: 
!  - the amg_c_slu_solver_type data structure containing the ingredients
!    to interface with the SuperLU package. 
!    1. The factorization is restricted to the diagonal block of the
!       current image.
!
module amg_c_slu_solver

  use iso_c_binding
  use amg_c_base_solver_mod

#if defined(IPK8)

  type, extends(amg_c_base_solver_type) :: amg_c_slu_solver_type

  end type amg_c_slu_solver_type

#else  

  type, extends(amg_c_base_solver_type) :: amg_c_slu_solver_type
    type(c_ptr)                 :: lufactors=c_null_ptr
    integer(c_long_long)        :: symbsize=0, numsize=0
  contains
    procedure, pass(sv) :: build   => c_slu_solver_bld
    procedure, pass(sv) :: apply_a => c_slu_solver_apply
    procedure, pass(sv) :: apply_v => c_slu_solver_apply_vect
    procedure, pass(sv) :: free    => c_slu_solver_free
    procedure, pass(sv) :: clear_data  => c_slu_solver_clear_data
    procedure, pass(sv) :: descr   => c_slu_solver_descr
    procedure, pass(sv) :: sizeof  => c_slu_solver_sizeof
    procedure, nopass   :: get_fmt => c_slu_solver_get_fmt
    procedure, nopass   :: get_id  => c_slu_solver_get_id
    final               :: c_slu_solver_finalize
  end type amg_c_slu_solver_type


  private :: c_slu_solver_bld, c_slu_solver_apply, &
       &  c_slu_solver_free,   c_slu_solver_descr, &
       &  c_slu_solver_sizeof, c_slu_solver_apply_vect, &
       &  c_slu_solver_get_fmt, c_slu_solver_get_id, &
       &  c_slu_solver_clear_data
  private :: c_slu_solver_finalize



  interface 
    function amg_cslu_fact(n,nnz,values,rowptr,colind,&
         & lufactors)&
         & bind(c,name='amg_cslu_fact') result(info)
      use iso_c_binding
      integer(c_int), value :: n,nnz
      integer(c_int)        :: info
      integer(c_int)        :: rowptr(*),colind(*)
      complex(c_float_complex)        :: values(*)
      type(c_ptr)           :: lufactors
    end function amg_cslu_fact
  end interface

  interface 
    function amg_cslu_solve(itrans,n,nrhs,b,ldb,lufactors)&
         & bind(c,name='amg_cslu_solve') result(info)
      use iso_c_binding
      integer(c_int)        :: info
      integer(c_int), value :: itrans,n,nrhs,ldb
      complex(c_float_complex)        :: b(ldb,*)
      type(c_ptr), value    :: lufactors
    end function amg_cslu_solve
  end interface

  interface 
    function amg_cslu_free(lufactors)&
         & bind(c,name='amg_cslu_free') result(info)
      use iso_c_binding
      integer(c_int)        :: info
      type(c_ptr), value    :: lufactors
    end function amg_cslu_free
  end interface

contains

  subroutine c_slu_solver_apply(alpha,sv,x,beta,y,desc_data,&
       & trans,work,info,init,initu)
    use psb_base_mod
    implicit none 
    type(psb_desc_type), intent(in)      :: desc_data
    class(amg_c_slu_solver_type), intent(inout) :: sv
    complex(psb_spk_),intent(inout)         :: x(:)
    complex(psb_spk_),intent(inout)         :: y(:)
    complex(psb_spk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    complex(psb_spk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info
    character, intent(in), optional       :: init
    complex(psb_spk_),intent(inout), optional :: initu(:)

    integer    :: n_row,n_col
    complex(psb_spk_), pointer :: ww(:)
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='c_slu_solver_apply'

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
             & a_err='complex(psb_spk_)')
        goto 9999      
      end if
    endif

    ww(1:n_row) = x(1:n_row)
    select case(trans_)
    case('N')
      info = amg_cslu_solve(0,n_row,1,ww,n_row,sv%lufactors)
    case('T')
      info = amg_cslu_solve(1,n_row,1,ww,n_row,sv%lufactors)
    case('C')
      info = amg_cslu_solve(2,n_row,1,ww,n_row,sv%lufactors)
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

  end subroutine c_slu_solver_apply
  
  subroutine c_slu_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
       & trans,work,wv,info,init,initu)
    use psb_base_mod
    implicit none 
    type(psb_desc_type), intent(in)      :: desc_data
    class(amg_c_slu_solver_type), intent(inout) :: sv
    type(psb_c_vect_type),intent(inout)  :: x
    type(psb_c_vect_type),intent(inout)  :: y
    complex(psb_spk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)           :: trans
    complex(psb_spk_),target, intent(inout) :: work(:)
    type(psb_c_vect_type),intent(inout) :: wv(:)
    integer, intent(out)                  :: info
    character, intent(in), optional                :: init
    type(psb_c_vect_type),intent(inout), optional   :: initu

    integer    :: err_act
    character(len=20)  :: name='c_slu_solver_apply_vect'

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
  
  end subroutine c_slu_solver_apply_vect

  subroutine c_slu_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

    use psb_base_mod

    Implicit None

    ! Arguments
    type(psb_cspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(inout)                  :: desc_a 
    class(amg_c_slu_solver_type), intent(inout)         :: sv
    integer, intent(out)                                :: info
    type(psb_cspmat_type), intent(in), target, optional :: b
    class(psb_c_base_sparse_mat), intent(in), optional  :: amold
    class(psb_c_base_vect_type), intent(in), optional   :: vmold
    class(psb_i_base_vect_type), intent(in), optional  :: imold
    ! Local variables
    type(psb_cspmat_type) :: atmp
    type(psb_c_csc_sparse_mat) :: acsc
    type(psb_c_coo_sparse_mat) :: acoo
    integer :: n_row,n_col, nrow_a, nztota
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='c_slu_solver_bld', ch_err
    
    info=psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = desc_a%get_context()
    call psb_info(ictxt, me, np)
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
    info = amg_cslu_fact(nrow_a,nztota,acsc%val,&
         & acsc%icp,acsc%ia,sv%lufactors)
    
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='amg_cslu_fact'
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
  end subroutine c_slu_solver_bld

  subroutine c_slu_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(amg_c_slu_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='c_slu_solver_free'

    call psb_erractionsave(err_act)

    info = psb_success_ 

    call sv%clear_data(info)

    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
  return
  end subroutine c_slu_solver_free

  subroutine c_slu_solver_clear_data(sv,info)

    Implicit None

    ! Arguments
    class(amg_c_slu_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='c_slu_solver_clear_data'

    call psb_erractionsave(err_act)

    info = psb_success_ 
    if (c_associated(sv%lufactors)) info = amg_cslu_free(sv%lufactors)
    sv%lufactors = c_null_ptr
    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine c_slu_solver_clear_data

  subroutine c_slu_solver_finalize(sv)

    Implicit None

    ! Arguments
    type(amg_c_slu_solver_type), intent(inout) :: sv
    integer :: info
    Integer :: err_act
    character(len=20)  :: name='c_slu_solver_finalize'

    call sv%free(info) 

    return
  
  end subroutine c_slu_solver_finalize

  subroutine c_slu_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(amg_c_slu_solver_type), intent(in) :: sv
    integer, intent(out)                     :: info
    integer, intent(in), optional            :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='amg_c_slu_solver_descr'
    integer :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = psb_out_unit
    endif
    
    write(iout_,*) '  SuperLU Sparse Factorization Solver. '

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine c_slu_solver_descr

  function c_slu_solver_sizeof(sv) result(val)

    implicit none 
    ! Arguments
    class(amg_c_slu_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer             :: i

    val = 2*psb_sizeof_ip + psb_sizeof_dp
    val = val + sv%symbsize
    val = val + sv%numsize
    return
  end function c_slu_solver_sizeof

  function c_slu_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "SuperLU solver"
  end function c_slu_solver_get_fmt

  function c_slu_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_slu_
  end function c_slu_solver_get_id
#endif
end module amg_c_slu_solver

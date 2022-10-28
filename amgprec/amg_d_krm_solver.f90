!  
!   
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.7)
!    
!    (C) Copyright 2021 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Fabio Durastante 
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
! moved here from amg4psblas-extension
!
!
!                             AMG4PSBLAS  Extensions
!
!    (C) Copyright 2019
!
!                        Salvatore Filippone  Cranfield University
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
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
! File: amg_d_krm_solver_mod.f90
!
! Module: amg_d_krm_solver_mod
!
module amg_d_krm_solver

  use amg_d_base_solver_mod
  use amg_d_prec_type

  type, extends(amg_d_base_solver_type) :: amg_d_krm_solver_type
    !
    logical               :: global
    character(len=16)     :: method, kprec, sub_solve
    integer(psb_ipk_)     :: irst, istopc, itmax, itrace, i_sub_solve
    integer(psb_ipk_)     :: fillin
    real(psb_dpk_)      :: eps
    type(amg_dprec_type)  :: prec
    type(psb_desc_type)   :: desc_local
    type(psb_dspmat_type) :: a_local
    type(psb_d_vect_type) :: x_local, z_local
    type(psb_dspmat_type), pointer  :: a=>null()
  contains
    !
    !
    procedure, pass(sv) :: dump    => d_krm_solver_dmp
    procedure, pass(sv) :: check   => d_krm_solver_check
    procedure, pass(sv) :: clone   => d_krm_solver_clone
    procedure, pass(sv) :: clone_settings   => d_krm_solver_clone_settings
    procedure, pass(sv) :: cnv     => d_krm_solver_cnv
    procedure, pass(sv) :: apply_v => amg_d_krm_solver_apply_vect
    procedure, pass(sv) :: apply_a => amg_d_krm_solver_apply
    procedure, pass(sv) :: clear_data  => d_krm_solver_clear_data
    procedure, pass(sv) :: free    => d_krm_solver_free
    procedure, pass(sv) :: cseti   => d_krm_solver_cseti
    procedure, pass(sv) :: csetc   => d_krm_solver_csetc
    procedure, pass(sv) :: csetr   => d_krm_solver_csetr
    procedure, pass(sv) :: sizeof  => d_krm_solver_sizeof
    procedure, pass(sv) :: get_nzeros => d_krm_solver_get_nzeros
    !procedure, nopass   :: get_id  => d_krm_solver_get_id
    procedure, pass(sv) :: is_global => d_krm_solver_is_global
    procedure, nopass   :: is_iterative => d_krm_solver_is_iterative


    !
    ! These methods are specific for the new solver type
    ! and therefore need to be overridden
    !
    procedure, pass(sv) :: descr   => d_krm_solver_descr
    procedure, pass(sv) :: default => d_krm_solver_default
    procedure, pass(sv) :: build   => amg_d_krm_solver_bld
    procedure, nopass   :: get_fmt => d_krm_solver_get_fmt
  end type amg_d_krm_solver_type


  private ::  d_krm_solver_get_fmt, d_krm_solver_descr, d_krm_solver_default

  interface
    subroutine amg_d_krm_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, amg_d_krm_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type, psb_ipk_
      implicit none
      type(psb_desc_type), intent(in)             :: desc_data
      class(amg_d_krm_solver_type), intent(inout) :: sv
      type(psb_d_vect_type),intent(inout)         :: x
      type(psb_d_vect_type),intent(inout)         :: y
      real(psb_dpk_),intent(in)                    :: alpha,beta
      character(len=1),intent(in)                   :: trans
      real(psb_dpk_),target, intent(inout)         :: work(:)
      type(psb_d_vect_type),intent(inout)         :: wv(:)
      integer(psb_ipk_), intent(out)                :: info
      character, intent(in), optional                :: init
      type(psb_d_vect_type),intent(inout), optional   :: initu
    end subroutine amg_d_krm_solver_apply_vect
  end interface

  interface
    subroutine amg_d_krm_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, amg_d_krm_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type, psb_ipk_
      implicit none
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_d_krm_solver_type), intent(inout) :: sv
      real(psb_dpk_),intent(inout)         :: x(:)
      real(psb_dpk_),intent(inout)         :: y(:)
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      real(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_d_krm_solver_apply
  end interface

  interface
    subroutine amg_d_krm_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, amg_d_krm_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none
      type(psb_dspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a
      class(amg_d_krm_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_dspmat_type), intent(in), target, optional :: b
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_d_krm_solver_bld
  end interface


contains

  !
  !
  subroutine d_krm_solver_default(sv)

    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv

    sv%method     = 'bicgstab'
    sv%kprec      = 'bjac'
    sv%sub_solve  = 'ilu'
    sv%i_sub_solve  = -1
    sv%fillin = 1
    sv%irst   = 30
    sv%istopc = 2
    sv%itmax  = 40
    sv%itrace = -1
    sv%eps    = 1.d-6
    sv%global = .false.

    return
  end subroutine d_krm_solver_default

  function d_krm_solver_get_nzeros(sv) result(val)

    implicit none
    ! Arguments
    class(amg_d_krm_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val

    val = sv%prec%get_nzeros()

    return
  end function d_krm_solver_get_nzeros

  function d_krm_solver_sizeof(sv) result(val)

    implicit none
    ! Arguments
    class(amg_d_krm_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val

    val = sv%prec%sizeof() + sv%desc_local%sizeof() + sv%a_local%sizeof()

    return
  end function d_krm_solver_sizeof


  subroutine d_krm_solver_check(sv,info)

    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_krm_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_


    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine d_krm_solver_check

  subroutine d_krm_solver_cseti(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_krm_solver_cseti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(psb_toupper(trim(what)))
    case('KRM_IRST')
      sv%irst = val
    case('KRM_ISTOPC')
      sv%istopc = val
    case('KRM_ITMAX')
      sv%itmax = val
    case('KRM_ITRACE')
      sv%itrace = val
    case('KRM_SUB_SOLVE')
      sv%i_sub_solve = val
    case('KRM_FILLIN')
      sv%fillin = val
    case default
      call sv%amg_d_base_solver_type%set(what,val,info,idx=idx)
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_krm_solver_cseti

  subroutine d_krm_solver_csetc(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    character(len=*), intent(in)                  :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act, ival
    character(len=20)  :: name='d_krm_solver_csetc'

    info = psb_success_
    call psb_erractionsave(err_act)


    select case(psb_toupper(trim(what)))
    case('KRM_METHOD')
      sv%method = psb_toupper(trim(val))
    case('KRM_KPREC')
      sv%kprec = psb_toupper(trim(val))
    case('KRM_SUB_SOLVE')
      sv%sub_solve = psb_toupper(trim(val))
    case('KRM_GLOBAL')
      select case(psb_toupper(trim(val)))
      case('LOCAL','FALSE')
        sv%global = .false.
      case('GLOBAL','TRUE')
        sv%global =.true.
      end select
    case default
      call sv%amg_d_base_solver_type%set(what,val,info,idx=idx)
    end select


    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_krm_solver_csetc

  subroutine d_krm_solver_csetr(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what
    real(psb_dpk_), intent(in)                  :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_krm_solver_csetr'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case(psb_toupper(what))
    case('KRM_EPS')
      sv%eps = val
    case default
      call sv%amg_d_base_solver_type%set(what,val,info,idx=idx)
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_krm_solver_csetr

  subroutine d_krm_solver_clear_data(sv,info)
    use psb_base_mod, only : psb_exit
    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_) :: err_act
    type(psb_ctxt_type)  :: l_ctxt
    character(len=20)  :: name='d_krm_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    call sv%prec%free(info)
    call sv%a_local%free()
    call sv%x_local%free(info)
    call sv%z_local%free(info)
    if ((.not.sv%global).and.sv%desc_local%is_ok()) then
      l_ctxt=sv%desc_local%get_context()
      call psb_exit(l_ctxt,close=.false.)
    end if
    call sv%desc_local%free(info)
    nullify(sv%a)
    call psb_erractionrestore(err_act)
    return
  end subroutine d_krm_solver_clear_data


  subroutine d_krm_solver_free(sv,info)
    use psb_base_mod, only : psb_exit
    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    type(psb_ctxt_type) :: l_ctxt
    character(len=20)  :: name='d_krm_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    call sv%clear_data(info)

    call psb_erractionrestore(err_act)
    return
  end subroutine d_krm_solver_free

  function d_krm_solver_get_fmt() result(val)
    implicit none
    character(len=32)  :: val

    val = "KRM solver"
  end function d_krm_solver_get_fmt

  subroutine d_krm_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)             :: info
    integer(psb_ipk_), intent(in), optional    :: iout
    logical, intent(in), optional              :: coarse
    character(len=*), intent(in), optional  :: prefix

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='amg_d_krm_solver_descr'
    integer(psb_ipk_) :: iout_
    character(1024)    :: prefix_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then
      iout_ = iout
    else
      iout_ = psb_out_unit
    endif
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if

    if (sv%global) then
      write(iout_,*) trim(prefix_), '  Krylov solver (global)'
    else
      write(iout_,*) trim(prefix_), '  Krylov solver (local) '
    end if
    write(iout_,*) trim(prefix_), '    method: ',sv%method
    write(iout_,*) trim(prefix_), '     kprec: ',sv%kprec
    call sv%prec%descr(info,iout_,prefix='KRM : '//prefix_)
    write(iout_,*) trim(prefix_), '     itmax: ',sv%itmax
    write(iout_,*) trim(prefix_), '       eps: ',sv%eps

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_krm_solver_descr

  subroutine d_krm_solver_cnv(sv,info,amold,vmold,imold)
    implicit none
    class(amg_d_krm_solver_type), intent(inout)         :: sv
    integer(psb_ipk_), intent(out)                      :: info
    class(psb_d_base_sparse_mat), intent(in), optional  :: amold
    class(psb_d_base_vect_type), intent(in), optional   :: vmold
    class(psb_i_base_vect_type), intent(in), optional   :: imold

    call sv%prec%cnv(info,amold=amold,vmold=vmold,imold=imold)

  end subroutine d_krm_solver_cnv

  subroutine d_krm_solver_clone(sv,svout,info)
    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout)               :: sv
    class(amg_d_base_solver_type), allocatable, intent(inout) :: svout
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_

    call svout%free(info)
    allocate(svout,stat=info,mold=sv)
    select type(so=>svout)
    class is(amg_d_krm_solver_type)
      so%method = sv%method
      so%kprec  = sv%kprec
      so%sub_solve = sv%sub_solve
      so%i_sub_solve = sv%i_sub_solve
      so%fillin = sv%fillin
      so%irst   = sv%irst
      so%istopc = sv%istopc
      so%itmax  = sv%itmax
      so%itrace = sv%itrace
      so%global = sv%global
      so%eps    = sv%eps
      call sv%a_local%clone(so%a_local,info)
      call sv%desc_local%clone(so%desc_local,info)
      call sv%prec%clone(so%prec,info)
    class default
      info = psb_err_internal_error_
    end select

  end subroutine d_krm_solver_clone


  subroutine d_krm_solver_clone_settings(sv,svout,info)
    Implicit None

    ! Arguments
    class(amg_d_krm_solver_type), intent(inout)               :: sv
    class(amg_d_base_solver_type), intent(inout) :: svout
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_

    select type(so=>svout)
    class is(amg_d_krm_solver_type)
      so%method = sv%method
      so%kprec  = sv%kprec
      so%sub_solve = sv%sub_solve
      so%i_sub_solve = sv%i_sub_solve
      so%fillin = sv%fillin
      so%irst   = sv%irst
      so%istopc = sv%istopc
      so%itmax  = sv%itmax
      so%itrace = sv%itrace
      so%global = sv%global
      so%eps    = sv%eps
    class default
      info = psb_err_internal_error_
    end select

  end subroutine d_krm_solver_clone_settings

  subroutine d_krm_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)
    implicit none
    class(amg_d_krm_solver_type), intent(in) :: sv
    type(psb_desc_type), intent(in)             :: desc
    integer(psb_ipk_), intent(in)              :: level
    integer(psb_ipk_), intent(out)             :: info
    character(len=*), intent(in), optional     :: prefix, head
    logical, optional, intent(in)              :: solver, global_num


    call sv%prec%dump(info,prefix=prefix,head=head)

  end subroutine d_krm_solver_dmp
  !
  ! Notify whether KRM is used as a global solver
  !
  function d_krm_solver_is_global(sv) result(val)
    implicit none
    class(amg_d_krm_solver_type), intent(in) :: sv
    logical  :: val

    val = (sv%global)
  end function d_krm_solver_is_global
  !
  function d_krm_solver_is_iterative() result(val)
    implicit none
    logical  :: val

    val = .true.
  end function d_krm_solver_is_iterative

end module amg_d_krm_solver

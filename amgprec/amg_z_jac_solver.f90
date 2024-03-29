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
!
! File: amg_z_jac_solver_mod.f90
!
! Module: amg_z_jac_solver_mod
!
!  This module defines: 
!  - the amg_z_jac_solver_type data structure containing the ingredients
!    for a local Jacobi  iteration. The iterations are local to a process
!    (they operate on the block diagonal).
!    
!
module amg_z_jac_solver

  use amg_z_base_solver_mod

  type, extends(amg_z_base_solver_type) :: amg_z_jac_solver_type
    type(psb_zspmat_type)      :: a
    type(psb_z_vect_type), allocatable :: dv
    complex(psb_dpk_), allocatable        :: d(:)
    integer(psb_ipk_)          :: sweeps
    real(psb_dpk_)             :: eps
  contains
    procedure, pass(sv) :: dump    => amg_z_jac_solver_dmp
    procedure, pass(sv) :: check   => z_jac_solver_check
    procedure, pass(sv) :: clone   => amg_z_jac_solver_clone
    procedure, pass(sv) :: clone_settings => amg_z_jac_solver_clone_settings
    procedure, pass(sv) :: clear_data     => amg_z_jac_solver_clear_data
    procedure, pass(sv) :: build   => amg_z_jac_solver_bld
    procedure, pass(sv) :: cnv     => amg_z_jac_solver_cnv
    procedure, pass(sv) :: apply_v => amg_z_jac_solver_apply_vect
    procedure, pass(sv) :: apply_a => amg_z_jac_solver_apply
    procedure, pass(sv) :: free    => z_jac_solver_free
    procedure, pass(sv) :: cseti   => z_jac_solver_cseti
    procedure, pass(sv) :: csetc   => z_jac_solver_csetc
    procedure, pass(sv) :: csetr   => z_jac_solver_csetr
    procedure, pass(sv) :: descr   => z_jac_solver_descr
    procedure, pass(sv) :: default => z_jac_solver_default
    procedure, pass(sv) :: sizeof  => z_jac_solver_sizeof
    procedure, pass(sv) :: get_nzeros => z_jac_solver_get_nzeros
    procedure, nopass   :: get_wrksz  => z_jac_solver_get_wrksize
    procedure, nopass   :: get_fmt    => z_jac_solver_get_fmt
    procedure, nopass   :: get_id     => z_jac_solver_get_id
    procedure, nopass   :: is_iterative => z_jac_solver_is_iterative
  end type amg_z_jac_solver_type

  type, extends(amg_z_jac_solver_type) :: amg_z_l1_jac_solver_type
  contains
    procedure, pass(sv) :: build   => amg_z_l1_jac_solver_bld
    procedure, pass(sv) :: descr   => z_l1_jac_solver_descr
    procedure, nopass   :: get_fmt    => z_l1_jac_solver_get_fmt
    procedure, nopass   :: get_id     => z_l1_jac_solver_get_id
  end type amg_z_l1_jac_solver_type


  private :: z_jac_solver_bld, z_jac_solver_apply, &
       &  z_jac_solver_free,   &
       &  z_jac_solver_descr,  z_jac_solver_sizeof, &
       &  z_jac_solver_default, z_jac_solver_dmp, &
       &  z_jac_solver_apply_vect, z_jac_solver_get_nzeros, &
       &  z_jac_solver_get_fmt, z_jac_solver_check,&
       &  z_jac_solver_is_iterative, &
       &  z_jac_solver_get_id, z_jac_solver_get_wrksize

  interface 
    subroutine amg_z_jac_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, amg_z_jac_solver_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)             :: desc_data
      class(amg_z_jac_solver_type), intent(inout) :: sv
      type(psb_z_vect_type),intent(inout)         :: x
      type(psb_z_vect_type),intent(inout)         :: y
      complex(psb_dpk_),intent(in)                    :: alpha,beta
      character(len=1),intent(in)                   :: trans
      complex(psb_dpk_),target, intent(inout)         :: work(:)
      type(psb_z_vect_type),intent(inout)         :: wv(:)
      integer(psb_ipk_), intent(out)                :: info
      character, intent(in), optional                :: init
      type(psb_z_vect_type),intent(inout), optional   :: initu
    end subroutine amg_z_jac_solver_apply_vect
  end interface

  interface 
    subroutine amg_z_jac_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info,init,initu)
      import :: psb_desc_type, amg_z_jac_solver_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_z_jac_solver_type), intent(inout) :: sv
      complex(psb_dpk_),intent(inout)         :: x(:)
      complex(psb_dpk_),intent(inout)         :: y(:)
      complex(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      complex(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      complex(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_z_jac_solver_apply
  end interface

  interface 
    subroutine amg_z_jac_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, amg_z_jac_solver_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none 
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(amg_z_jac_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_zspmat_type), intent(in), target, optional :: b
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_z_jac_solver_bld
  end interface

  interface 
    subroutine amg_z_l1_jac_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, amg_z_l1_jac_solver_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none 
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(amg_z_l1_jac_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_zspmat_type), intent(in), target, optional :: b
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_z_l1_jac_solver_bld
  end interface

  interface 
    subroutine amg_z_jac_solver_cnv(sv,info,amold,vmold,imold)
      import :: amg_z_jac_solver_type, psb_dpk_, &
           & psb_z_base_sparse_mat, psb_z_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none 
      class(amg_z_jac_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_z_jac_solver_cnv
  end interface
  
  interface 
    subroutine amg_z_jac_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)
      import :: psb_desc_type, amg_z_jac_solver_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type, &
           & psb_ipk_
      implicit none 
      class(amg_z_jac_solver_type), intent(in) :: sv
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional     :: prefix, head
      logical, optional, intent(in)              :: solver, global_num
    end subroutine amg_z_jac_solver_dmp
  end interface
  
  interface
    subroutine amg_z_jac_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_solver_type, amg_z_jac_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_z_jac_solver_type), intent(inout)               :: sv
      class(amg_z_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)              :: info
    end subroutine amg_z_jac_solver_clone
  end interface
  
!!$  interface
!!$    subroutine amg_z_l1_jac_solver_clone(sv,svout,info)
!!$      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
!!$           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
!!$           & amg_z_base_solver_type, amg_z_l1_jac_solver_type, psb_ipk_
!!$      Implicit None
!!$      
!!$      ! Arguments
!!$      class(amg_z_l1_jac_solver_type), intent(inout)               :: sv
!!$      class(amg_z_base_solver_type), allocatable, intent(inout) :: svout
!!$      integer(psb_ipk_), intent(out)              :: info
!!$    end subroutine amg_z_l1_jac_solver_clone
!!$  end interface

  interface
    subroutine amg_z_jac_solver_clone_settings(sv,svout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_solver_type, amg_z_jac_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_z_jac_solver_type), intent(inout) :: sv
      class(amg_z_base_solver_type), intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_z_jac_solver_clone_settings
  end interface

  interface
    subroutine amg_z_jac_solver_clear_data(sv,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_jac_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_z_jac_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_z_jac_solver_clear_data
  end interface

contains

  subroutine z_jac_solver_default(sv)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(inout) :: sv

    sv%sweeps = ione
    sv%eps    = dzero
    
    return
  end subroutine z_jac_solver_default

  subroutine z_jac_solver_check(sv,info)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='z_jac_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call amg_check_def(sv%sweeps,&
         & 'Jacobi sweeps',ione,is_int_positive)

    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine z_jac_solver_check

  subroutine z_jac_solver_cseti(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(inout) :: sv 
    character(len=*), intent(in)                  :: what 
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='z_jac_solver_cseti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(psb_toupper(what))
    case('SOLVER_SWEEPS')
      sv%sweeps = val      
    case default
      call sv%amg_z_base_solver_type%set(what,val,info,idx=idx)
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_jac_solver_cseti

  subroutine z_jac_solver_csetc(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what 
    character(len=*), intent(in)                  :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act, ival
    character(len=20)  :: name='z_jac_solver_csetc'

    info = psb_success_
    call psb_erractionsave(err_act)


    call sv%amg_z_base_solver_type%set(what,val,info,idx=idx)
  
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_jac_solver_csetc
  
  subroutine z_jac_solver_csetr(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(inout) :: sv 
    character(len=*), intent(in)                  :: what 
    real(psb_dpk_), intent(in)                     :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='z_jac_solver_csetr'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case(psb_toupper(what))
    case('SOLVER_EPS') 
      sv%eps = val
    case default
      call sv%amg_z_base_solver_type%set(what,val,info,idx=idx)
    end select


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_jac_solver_csetr

  subroutine z_jac_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='z_jac_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    call sv%a%free()
    if (allocated(sv%dv)) then
      call sv%dv%free(info)
      deallocate(sv%dv)
    end if
    if (allocated(sv%d)) deallocate(sv%d)

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_jac_solver_free

  subroutine z_jac_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_z_jac_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: iout
    logical, intent(in), optional             :: coarse
    character(len=*), intent(in), optional    :: prefix

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='amg_z_jac_solver_descr'
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

    if (sv%eps<=dzero) then 
      write(iout_,*) trim(prefix_), '  Jacobi iterative solver with  ',&
           &  sv%sweeps,' sweeps'
    else
      write(iout_,*) trim(prefix_), '  Jacobi iterative solver with  tolerance',&
           &  sv%eps,' and maxit', sv%sweeps
    end if
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_jac_solver_descr

  function z_jac_solver_get_nzeros(sv) result(val)

    implicit none 
    ! Arguments
    class(amg_z_jac_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_)        :: i
    
    val = 0 
    val = val + sv%a%get_nzeros()
    val = val + sv%dv%get_nrows()

    return
  end function z_jac_solver_get_nzeros

  function z_jac_solver_sizeof(sv) result(val)

    implicit none 
    ! Arguments
    class(amg_z_jac_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_)        :: i

    val = psb_sizeof_ip 
    val = val + sv%a%sizeof()
    val = val + sv%dv%sizeof()

    return
  end function z_jac_solver_sizeof

  function z_jac_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Jacobi solver"
  end function z_jac_solver_get_fmt

  function z_jac_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_jac_
  end function z_jac_solver_get_id

  !
  ! If this is true, then the solver needs a starting
  ! guess. Currently only handled in JAC smoother. 
  ! 
  function z_jac_solver_is_iterative() result(val)
    implicit none 
    logical  :: val

    val = .true.
  end function z_jac_solver_is_iterative
    
  function z_jac_solver_get_wrksize() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = 2
  end function z_jac_solver_get_wrksize

  subroutine z_l1_jac_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_z_l1_jac_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: iout
    logical, intent(in), optional             :: coarse
    character(len=*), intent(in), optional    :: prefix

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='amg_z_l1_jac_solver_descr'
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

    if (sv%eps<=dzero) then 
      write(iout_,*) trim(prefix_), '  L1-Jacobi iterative solver with  ',&
           &  sv%sweeps,' sweeps'
    else
      write(iout_,*) trim(prefix_), '  L1-Jacobi iterative solver with  tolerance',&
           &  sv%eps,' and maxit', sv%sweeps
    end if
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_l1_jac_solver_descr

  function z_l1_jac_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "L1-Jacobi solver"
  end function z_l1_jac_solver_get_fmt

  function z_l1_jac_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_l1_jac_
  end function z_l1_jac_solver_get_id
  
end module amg_z_jac_solver

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
!        Daniela di Serafino
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
!         software without specific prior written permission.
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
! File: amg_d_richards_smoother.f90
!
! Module: amg_d_richards_smoother
!
!  This module defines:
!    the amg_d_richards_smoother_type data structure containing the
!    smoother for a preconditioned Richards iteration.
!  The smoother applies the iterative method:
!    x_{k+1} = x_k + omega * P^{-1} * (b - A * x_k)
!  where P is a preconditioner (the solver component) acting on the global matrix A.
!  This allows using distributed solvers like MUMPS on the full system.
!
module amg_d_richards_smoother

  use amg_d_base_smoother_mod

  type, extends(amg_d_base_smoother_type) :: amg_d_richards_smoother_type
    ! The local solver component is inherited from the
    ! parent type, but acts on the global matrix.
    !    class(amg_d_base_solver_type), allocatable :: sv
    !
    type(psb_dspmat_type), pointer  :: pa => null()
    integer(psb_lpk_)               :: global_nnz_tot
    logical                         :: checkres
    logical                         :: printres
    integer(psb_ipk_)               :: checkiter
    integer(psb_ipk_)               :: printiter
    real(psb_dpk_)                  :: tol
    real(psb_dpk_)                  :: omega
  contains
    procedure, pass(sm) :: apply_v => amg_d_richards_smoother_apply_vect
    procedure, pass(sm) :: apply_a => amg_d_richards_smoother_apply
    procedure, pass(sm) :: dump    => amg_d_richards_smoother_dmp
    procedure, pass(sm) :: build   => amg_d_richards_smoother_bld
    procedure, pass(sm) :: cnv     => amg_d_richards_smoother_cnv
    procedure, pass(sm) :: clone   => amg_d_richards_smoother_clone
    procedure, pass(sm) :: clone_settings => amg_d_richards_smoother_clone_settings
    procedure, pass(sm) :: clear_data     => amg_d_richards_smoother_clear_data
    procedure, pass(sm) :: free    => d_richards_smoother_free
    procedure, pass(sm) :: cseti   => amg_d_richards_smoother_cseti
    procedure, pass(sm) :: csetc   => amg_d_richards_smoother_csetc
    procedure, pass(sm) :: csetr   => amg_d_richards_smoother_csetr
    procedure, pass(sm) :: descr   => amg_d_richards_smoother_descr
    procedure, pass(sm) :: sizeof  => d_richards_smoother_sizeof
    procedure, pass(sm) :: default => d_richards_smoother_default
    procedure, pass(sm) :: get_nzeros => d_richards_smoother_get_nzeros
    procedure, pass(sm) :: get_wrksz => d_richards_smoother_get_wrksize
    procedure, nopass   :: get_fmt    => d_richards_smoother_get_fmt
    procedure, nopass   :: get_id     => d_richards_smoother_get_id
  end type amg_d_richards_smoother_type

  private :: d_richards_smoother_free, &
       & d_richards_smoother_sizeof, d_richards_smoother_get_nzeros, &
       & d_richards_smoother_get_fmt, d_richards_smoother_get_id, &
       & d_richards_smoother_get_wrksize

  interface
    subroutine amg_d_richards_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,&
         & sweeps,work,wv,info,init,initu)
      import :: psb_desc_type, amg_d_richards_smoother_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_

      type(psb_desc_type), intent(in)                 :: desc_data
      class(amg_d_richards_smoother_type), intent(inout) :: sm
      type(psb_d_vect_type),intent(inout)           :: x
      type(psb_d_vect_type),intent(inout)           :: y
      real(psb_dpk_),intent(in)                      :: alpha,beta
      character(len=1),intent(in)                     :: trans
      integer(psb_ipk_), intent(in)                   :: sweeps
      real(psb_dpk_),target, intent(inout)           :: work(:)
      type(psb_d_vect_type),intent(inout)           :: wv(:)
      integer(psb_ipk_), intent(out)                  :: info
      character, intent(in), optional                :: init
      type(psb_d_vect_type),intent(inout), optional   :: initu
    end subroutine amg_d_richards_smoother_apply_vect
  end interface

  interface
    subroutine amg_d_richards_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,&
         & sweeps,work,info,init,initu)
      import :: psb_desc_type, amg_d_richards_smoother_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type, &
           & psb_ipk_
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_d_richards_smoother_type), intent(inout) :: sm
      real(psb_dpk_),intent(inout)         :: x(:)
      real(psb_dpk_),intent(inout)         :: y(:)
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      integer(psb_ipk_), intent(in)         :: sweeps
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      real(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_d_richards_smoother_apply
  end interface

  interface
    subroutine amg_d_richards_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)
      import :: psb_desc_type, amg_d_richards_smoother_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      type(psb_dspmat_type), intent(inout), target        :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a
      class(amg_d_richards_smoother_type), intent(inout)       :: sm
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_d_richards_smoother_bld
  end interface

  interface
    subroutine amg_d_richards_smoother_cnv(sm,info,amold,vmold,imold)
      import :: amg_d_richards_smoother_type, psb_dpk_, &
           & psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      class(amg_d_richards_smoother_type), intent(inout)       :: sm
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_d_richards_smoother_cnv
  end interface

  interface
    subroutine amg_d_richards_smoother_dmp(sm,desc,level,info,prefix,head,smoother,solver,global_num)
      import :: psb_dspmat_type, psb_d_vect_type, psb_d_base_vect_type, &
           & psb_dpk_, amg_d_richards_smoother_type, psb_epk_, psb_desc_type, &
           & psb_ipk_
      implicit none
      class(amg_d_richards_smoother_type), intent(in) :: sm
      type(psb_desc_type), intent(in)               :: desc
      integer(psb_ipk_), intent(in)               :: level
      integer(psb_ipk_), intent(out)              :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: smoother, solver, global_num
    end subroutine amg_d_richards_smoother_dmp
  end interface

  interface
    subroutine amg_d_richards_smoother_clone(sm,smout,info)
      import :: amg_d_richards_smoother_type, psb_dpk_, &
           & amg_d_base_smoother_type, psb_ipk_
      class(amg_d_richards_smoother_type), intent(inout)               :: sm
      class(amg_d_base_smoother_type), allocatable, intent(inout) :: smout
      integer(psb_ipk_), intent(out)                :: info
    end subroutine amg_d_richards_smoother_clone
  end interface

  interface
    subroutine amg_d_richards_smoother_clone_settings(sm,smout,info)
      import :: amg_d_richards_smoother_type, psb_dpk_, &
           & amg_d_base_smoother_type, psb_ipk_
      class(amg_d_richards_smoother_type), intent(inout)  :: sm
      class(amg_d_base_smoother_type), intent(inout) :: smout
      integer(psb_ipk_), intent(out)                :: info
    end subroutine amg_d_richards_smoother_clone_settings
  end interface

  interface
    subroutine amg_d_richards_smoother_clear_data(sm,info)
      import :: amg_d_richards_smoother_type, psb_dpk_, &
           & amg_d_base_smoother_type, psb_ipk_
      class(amg_d_richards_smoother_type), intent(inout)               :: sm
      integer(psb_ipk_), intent(out)                :: info
    end subroutine amg_d_richards_smoother_clear_data
  end interface

  interface
    subroutine amg_d_richards_smoother_descr(sm,info,iout,coarse,prefix)
      import :: amg_d_richards_smoother_type, psb_ipk_
      class(amg_d_richards_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: iout
      logical, intent(in), optional                :: coarse
      character(len=*), intent(in), optional       :: prefix
    end subroutine amg_d_richards_smoother_descr
  end interface

  interface
    subroutine amg_d_richards_smoother_cseti(sm,what,val,info,idx)
      import :: psb_dspmat_type, psb_d_vect_type, psb_d_base_vect_type, &
           & psb_dpk_, amg_d_richards_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none
      class(amg_d_richards_smoother_type), intent(inout) :: sm
      character(len=*), intent(in)                   :: what
      integer(psb_ipk_), intent(in)                  :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine amg_d_richards_smoother_cseti
  end interface

  interface
    subroutine amg_d_richards_smoother_csetc(sm,what,val,info,idx)
      import :: psb_dspmat_type, psb_d_vect_type, psb_d_base_vect_type, &
           & psb_dpk_, amg_d_richards_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none
      class(amg_d_richards_smoother_type), intent(inout) :: sm
      character(len=*), intent(in)                   :: what
      character(len=*), intent(in)                   :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine amg_d_richards_smoother_csetc
  end interface

  interface
    subroutine amg_d_richards_smoother_csetr(sm,what,val,info,idx)
      import :: psb_dspmat_type, psb_d_vect_type, psb_d_base_vect_type, &
           & psb_dpk_, amg_d_richards_smoother_type, psb_epk_, psb_desc_type, psb_ipk_
      implicit none
      class(amg_d_richards_smoother_type), intent(inout) :: sm
      character(len=*), intent(in)                   :: what
      real(psb_dpk_), intent(in)                   :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine amg_d_richards_smoother_csetr
  end interface

contains

  subroutine d_richards_smoother_free(sm,info)
    Implicit None

    ! Arguments
    class(amg_d_richards_smoother_type), intent(inout) :: sm
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='d_richards_smoother_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sm%sv)) then
      call sm%sv%free(info)
      if (info == psb_success_) deallocate(sm%sv,stat=info)
      if (info /= psb_success_) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999
      end if
    end if
    sm%pa => null()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_richards_smoother_free

  function d_richards_smoother_sizeof(sm) result(val)
    implicit none
    ! Arguments
    class(amg_d_richards_smoother_type), intent(in) :: sm
    integer(psb_epk_) :: val
    integer(psb_ipk_)        :: i

    val = psb_sizeof_lp
    if (allocated(sm%sv)) val = val + sm%sv%sizeof()

    return
  end function d_richards_smoother_sizeof

  subroutine d_richards_smoother_default(sm)
    Implicit None

    ! Arguments
    class(amg_d_richards_smoother_type), intent(inout) :: sm

    !
    ! Default: Richards iteration with omega=1.0 and no residual check
    !
    sm%checkres = .false.
    sm%printres = .false.
    sm%checkiter = -1
    sm%printiter = -1
    sm%tol = 0
    sm%omega = 1.0d0

    if (allocated(sm%sv)) then
      call sm%sv%default()
    end if

    return
  end subroutine d_richards_smoother_default

  function d_richards_smoother_get_nzeros(sm) result(val)
    implicit none
    ! Arguments
    class(amg_d_richards_smoother_type), intent(in) :: sm
    integer(psb_epk_) :: val
    integer(psb_ipk_)        :: i

    val = 0
    if (allocated(sm%sv)) val = val + sm%sv%get_nzeros()

    return
  end function d_richards_smoother_get_nzeros

  function d_richards_smoother_get_wrksize(sm) result(val)
    implicit none
    class(amg_d_richards_smoother_type), intent(inout) :: sm
    integer(psb_ipk_)  :: val

    val = 2
    if (allocated(sm%sv)) val = val + sm%sv%get_wrksz()

  end function d_richards_smoother_get_wrksize

  function d_richards_smoother_get_fmt() result(val)
    implicit none
    character(len=32)  :: val

    val = "Richards smoother"
  end function d_richards_smoother_get_fmt

  function d_richards_smoother_get_id() result(val)
    implicit none
    integer(psb_ipk_)  :: val

    val = amg_richardson_
  end function d_richards_smoother_get_id

end module amg_d_richards_smoother

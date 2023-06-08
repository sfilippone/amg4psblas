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
! File: amg_s_diag_solver_mod.f90
!
! Module: amg_s_diag_solver_mod
!
!  This module defines: 
!  - the amg_s_diag_solver_type data structure containing the 
!    simple diagonal solver. This extracts the main diagonal of a matrix
!    and precomputes its inverse. Combined with a Jacobi "smoother" generates
!    what are commonly known as the classic Jacobi iterations
!
module amg_s_diag_solver

  use amg_s_base_solver_mod

  type, extends(amg_s_base_solver_type) :: amg_s_diag_solver_type
    type(psb_s_vect_type), allocatable :: dv
    real(psb_spk_), allocatable        :: d(:)
  contains
    procedure, pass(sv) :: dump    => amg_s_diag_solver_dmp
    procedure, pass(sv) :: build   => amg_s_diag_solver_bld
    procedure, pass(sv) :: cnv     => amg_s_diag_solver_cnv
    procedure, pass(sv) :: clone   => amg_s_diag_solver_clone
    procedure, pass(sv) :: clear_data  => amg_s_diag_solver_clear_data
    procedure, pass(sv) :: apply_v => amg_s_diag_solver_apply_vect
    procedure, pass(sv) :: apply_a => amg_s_diag_solver_apply
    procedure, pass(sv) :: free    => s_diag_solver_free
    procedure, pass(sv) :: descr   => s_diag_solver_descr
    procedure, pass(sv) :: sizeof  => s_diag_solver_sizeof
    procedure, pass(sv) :: get_nzeros  => s_diag_solver_get_nzeros
    procedure, nopass   :: get_fmt   => s_diag_solver_get_fmt
    procedure, nopass   :: get_id    => s_diag_solver_get_id
  end type amg_s_diag_solver_type


  private :: s_diag_solver_free,  s_diag_solver_descr, &
       & s_diag_solver_sizeof, s_diag_solver_get_nzeros, &
       & s_diag_solver_get_fmt, s_diag_solver_get_id


  interface 
    subroutine amg_s_diag_solver_apply_vect(alpha,sv,x,beta,y,desc_data,& 
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
       & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
       & amg_s_diag_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)                :: desc_data
      class(amg_s_diag_solver_type), intent(inout) :: sv
      type(psb_s_vect_type), intent(inout)         :: x
      type(psb_s_vect_type), intent(inout)         :: y
      real(psb_spk_),intent(in)                     :: alpha,beta
      character(len=1),intent(in)                    :: trans
      real(psb_spk_),target, intent(inout)          :: work(:)
      type(psb_s_vect_type),intent(inout)          :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_s_vect_type),intent(inout), optional   :: initu
    end subroutine amg_s_diag_solver_apply_vect
  end interface
  
  interface 
    subroutine amg_s_diag_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
       & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
       & amg_s_diag_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)            :: desc_data
      class(amg_s_diag_solver_type), intent(inout) :: sv
      real(psb_spk_), intent(inout)             :: x(:)
      real(psb_spk_), intent(inout)             :: y(:)
      real(psb_spk_),intent(in)                 :: alpha,beta
      character(len=1),intent(in)                :: trans
      real(psb_spk_),target, intent(inout)      :: work(:)
      integer(psb_ipk_), intent(out)             :: info
      character, intent(in), optional       :: init
      real(psb_spk_),intent(inout), optional :: initu(:)
    end subroutine amg_s_diag_solver_apply
  end interface
  
  interface 
    subroutine amg_s_diag_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_diag_solver_type, psb_ipk_, psb_i_base_vect_type      
      type(psb_sspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                    :: desc_a 
      class(amg_s_diag_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                        :: info
      type(psb_sspmat_type), intent(in), target, optional :: b
      class(psb_s_base_sparse_mat), intent(in), optional  :: amold
      class(psb_s_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_s_diag_solver_bld
  end interface
  
  interface 
    subroutine amg_s_diag_solver_cnv(sv,info,amold,vmold,imold)
      import :: psb_s_base_sparse_mat, psb_s_base_vect_type, psb_spk_, &
           & amg_s_diag_solver_type, psb_ipk_, psb_i_base_vect_type      
      class(amg_s_diag_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                        :: info
      class(psb_s_base_sparse_mat), intent(in), optional  :: amold
      class(psb_s_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_s_diag_solver_cnv
  end interface

  interface 
    subroutine amg_s_diag_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)
      import :: psb_desc_type, amg_s_diag_solver_type, psb_s_vect_type, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type, &
           & psb_ipk_
      implicit none 
      class(amg_s_diag_solver_type), intent(in) :: sv
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional     :: prefix, head
      logical, optional, intent(in)              :: solver, global_num
    end subroutine amg_s_diag_solver_dmp
  end interface
  
  interface
    subroutine amg_s_diag_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, amg_s_diag_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_diag_solver_type), intent(inout)              :: sv
      class(amg_s_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_s_diag_solver_clone
  end interface

  interface
    subroutine amg_s_diag_solver_clear_data(sv,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_diag_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_diag_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_s_diag_solver_clear_data
  end interface
  
  
contains

  subroutine s_diag_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(amg_s_diag_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                 :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_diag_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sv%dv)) call sv%dv%free(info)
    
    if (allocated(sv%d)) then 
      deallocate(sv%d,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine s_diag_solver_free

  subroutine s_diag_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_s_diag_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), intent(in), optional     :: iout
    logical, intent(in), optional               :: coarse
    character(len=*), intent(in), optional      :: prefix

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='amg_s_diag_solver_descr'
    integer(psb_ipk_) :: iout_
    character(1024)    :: prefix_

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

    write(iout_,*) trim(prefix_), '  Diagonal local solver '

    return

  end subroutine s_diag_solver_descr

  function s_diag_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(amg_s_diag_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i

    val = 0
    if (allocated(sv%dv)) val = val + sv%dv%sizeof()

    return
  end function s_diag_solver_sizeof

  function s_diag_solver_get_nzeros(sv) result(val)
    implicit none 
    ! Arguments
    class(amg_s_diag_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i

    val = 0
    if (allocated(sv%dv)) val = val +  sv%dv%get_nrows()

    return
  end function s_diag_solver_get_nzeros

  function s_diag_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Diag solver"
  end function s_diag_solver_get_fmt

  function s_diag_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_diag_scale_
  end function s_diag_solver_get_id

end module amg_s_diag_solver

!
! Module: amg_s_l1_diag_solver_mod
!
!  This module defines: 
!  - the amg_s_l1_diag_solver_type data structure containing the 
!    L1  diagonal solver. 
!    The solver is defined as a diagonal containing in each element the
!    inverse of the sum of the absolute values of the matrix entries
!    along the corresponding row. 
!    Combined with a Jacobi "smoother" generates 
!    what are commonly known as the L1-Jacobi iterations
!

module amg_s_l1_diag_solver

  use amg_s_diag_solver

  type, extends(amg_s_diag_solver_type) :: amg_s_l1_diag_solver_type
  contains
    procedure, pass(sv) :: dump    => amg_s_l1_diag_solver_dmp
    procedure, pass(sv) :: build   => amg_s_l1_diag_solver_bld
    procedure, pass(sv) :: descr   => s_l1_diag_solver_descr
    procedure, nopass   :: get_fmt   => s_l1_diag_solver_get_fmt
    procedure, nopass   :: get_id    => s_l1_diag_solver_get_id
  end type amg_s_l1_diag_solver_type


  private :: s_l1_diag_solver_descr, &
       & s_l1_diag_solver_get_fmt, s_l1_diag_solver_get_id

  interface 
    subroutine amg_s_l1_diag_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_l1_diag_solver_type, psb_ipk_, psb_i_base_vect_type      
      type(psb_sspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                    :: desc_a 
      class(amg_s_l1_diag_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                        :: info
      type(psb_sspmat_type), intent(in), target, optional :: b
      class(psb_s_base_sparse_mat), intent(in), optional  :: amold
      class(psb_s_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_s_l1_diag_solver_bld
  end interface
  
  interface 
    subroutine amg_s_l1_diag_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)
      import :: psb_desc_type, amg_s_l1_diag_solver_type, psb_s_vect_type, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type, &
           & psb_ipk_
      implicit none 
      class(amg_s_l1_diag_solver_type), intent(in) :: sv
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional     :: prefix, head
      logical, optional, intent(in)              :: solver, global_num
    end subroutine amg_s_l1_diag_solver_dmp
  end interface
  
contains

  subroutine s_l1_diag_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_s_l1_diag_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), intent(in), optional     :: iout
    logical, intent(in), optional               :: coarse
    character(len=*), intent(in), optional  :: prefix

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='amg_s_l1_diag_solver_descr'
    integer(psb_ipk_) :: iout_
    character(1024)    :: prefix_

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

    write(iout_,*) trim(prefix_), '  L1 Diagonal solver '

    return

  end subroutine s_l1_diag_solver_descr

  function s_l1_diag_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "L1 Diag solver"
  end function s_l1_diag_solver_get_fmt

  function s_l1_diag_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_l1_diag_scale_
  end function s_l1_diag_solver_get_id

end module amg_s_l1_diag_solver


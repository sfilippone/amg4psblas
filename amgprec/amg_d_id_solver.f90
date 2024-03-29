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
!
! Identity solver. Reference for nullprec. 
!
!
module amg_d_id_solver

  use amg_d_base_solver_mod

  type, extends(amg_d_base_solver_type) :: amg_d_id_solver_type
  contains
    procedure, pass(sv) :: build   => d_id_solver_bld
    procedure, pass(sv) :: clone   => amg_d_id_solver_clone
    procedure, pass(sv) :: apply_v => amg_d_id_solver_apply_vect
    procedure, pass(sv) :: apply_a => amg_d_id_solver_apply
    procedure, pass(sv) :: free    => d_id_solver_free
    procedure, pass(sv) :: descr   => d_id_solver_descr
    procedure, nopass   :: get_fmt   => d_id_solver_get_fmt
    procedure, nopass   :: get_id    => d_id_solver_get_id
  end type amg_d_id_solver_type


  private :: d_id_solver_bld, &
       &  d_id_solver_free, d_id_solver_get_fmt, &
       &  d_id_solver_descr, d_id_solver_get_id

  interface 
    subroutine amg_d_id_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, & 
           & amg_d_id_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)              :: desc_data
      class(amg_d_id_solver_type), intent(inout) :: sv
      type(psb_d_vect_type),intent(inout)        :: x
      type(psb_d_vect_type),intent(inout)        :: y
      real(psb_dpk_),intent(in)                   :: alpha,beta
      character(len=1),intent(in)                  :: trans
      real(psb_dpk_),target, intent(inout)        :: work(:)
      type(psb_d_vect_type),intent(inout)        :: wv(:)
      integer(psb_ipk_), intent(out)               :: info
      character, intent(in), optional                :: init
      type(psb_d_vect_type),intent(inout), optional   :: initu
    end subroutine amg_d_id_solver_apply_vect
  end interface
  
  interface 
    subroutine amg_d_id_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, &
           & amg_d_id_solver_type, psb_ipk_
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_d_id_solver_type), intent(inout) :: sv
      real(psb_dpk_),intent(inout)         :: x(:)
      real(psb_dpk_),intent(inout)         :: y(:)
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional       :: init
      real(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_d_id_solver_apply
  end interface

  interface
    subroutine amg_d_id_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, &
           & amg_d_base_solver_type, amg_d_id_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_d_id_solver_type), intent(inout)                :: sv
      class(amg_d_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)             :: info
    end subroutine amg_d_id_solver_clone
  end interface

contains


  subroutine d_id_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

    Implicit None

    ! Arguments
    type(psb_dspmat_type), intent(in), target           :: a
    Type(psb_desc_type), Intent(inout)                  :: desc_a 
    class(amg_d_id_solver_type), intent(inout)          :: sv
    integer(psb_ipk_), intent(out)                      :: info
    type(psb_dspmat_type), intent(in), target, optional :: b
    class(psb_d_base_sparse_mat), intent(in), optional  :: amold
    class(psb_d_base_vect_type), intent(in), optional   :: vmold
    class(psb_i_base_vect_type), intent(in), optional   :: imold
    ! Local variables
    integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota
    real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer(psb_ipk_) :: i, err_act, debug_unit, debug_level
    character(len=20)  :: name='d_id_solver_bld', ch_err
    
    info=psb_success_

    return
  end subroutine d_id_solver_bld

  subroutine d_id_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(amg_d_id_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)               :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_id_solver_free'

    info = psb_success_

    return
  end subroutine d_id_solver_free

  subroutine d_id_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_d_id_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: iout
    logical, intent(in), optional             :: coarse
    character(len=*), intent(in), optional    :: prefix

    ! Local variables
    integer(psb_ipk_)      :: err_act
    character(len=20), parameter :: name='amg_d_id_solver_descr'
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
    
    write(iout_,*) trim(prefix_), '  Identity local solver '

    return

  end subroutine d_id_solver_descr

  function d_id_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Identity solver"
  end function d_id_solver_get_fmt

  function d_id_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_f_none_
  end function d_id_solver_get_id

end module amg_d_id_solver

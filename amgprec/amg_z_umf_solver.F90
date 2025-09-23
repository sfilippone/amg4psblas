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
! File: amg_z_umf_solver_mod.f90
!
! Module: amg_z_umf_solver_mod
!
!  This module defines: 
!  - the amg_z_umf_solver_type data structure containing the ingredients
!    to interface with the UMFPACK package. 
!    1. The factorization is restricted to the diagonal block of the
!       current image.
!
module amg_z_umf_solver

  use iso_c_binding
  use amg_z_base_solver_mod

#if defined(PSB_IPK8)
  type, extends(amg_z_base_solver_type) :: amg_z_umf_solver_type

  end type amg_z_umf_solver_type

#else
  
  type, extends(amg_z_base_solver_type) :: amg_z_umf_solver_type
    type(c_ptr)                 :: symbolic=c_null_ptr, numeric=c_null_ptr
    integer(c_long_long)        :: symbsize=0, numsize=0
  contains
    procedure, pass(sv) :: build   => amg_z_umf_solver_bld
    procedure, pass(sv) :: apply_a => amg_z_umf_solver_apply
    procedure, pass(sv) :: apply_v => amg_z_umf_solver_apply_vect
    procedure, pass(sv) :: free    => z_umf_solver_free
    procedure, pass(sv) :: clear_data  => z_umf_solver_clear_data
    procedure, pass(sv) :: descr   => z_umf_solver_descr
    procedure, pass(sv) :: sizeof  => z_umf_solver_sizeof
    procedure, nopass   :: get_fmt => z_umf_solver_get_fmt
    procedure, nopass   :: get_id  => z_umf_solver_get_id
    final               :: z_umf_solver_finalize
  end type amg_z_umf_solver_type


  private :: z_umf_solver_free,   z_umf_solver_descr, &
       &  z_umf_solver_sizeof, &
       &  z_umf_solver_get_fmt, z_umf_solver_get_id, &
       &  z_umf_solver_clear_data
  private :: z_umf_solver_finalize



  interface 
    function amg_zumf_fact(n,nnz,values,rowind,colptr,&
         & symptr,numptr,ssize,nsize)&
         & bind(c,name='amg_zumf_fact') result(info)
      use iso_c_binding
      integer(c_int), value :: n,nnz
      integer(c_int)        :: info
      integer(c_long_long)  :: ssize, nsize
      integer(c_int)        :: rowind(*),colptr(*)
      complex(c_double_complex)  :: values(*)
      type(c_ptr)           :: symptr, numptr
    end function amg_zumf_fact
  end interface

  interface 
    function amg_zumf_solve(itrans,n,x, b, ldb, numptr)&
         & bind(c,name='amg_zumf_solve') result(info)
      use iso_c_binding
      integer(c_int)        :: info
      integer(c_int), value :: itrans,n,ldb
      complex(c_double_complex) :: x(*), b(ldb,*)
      type(c_ptr), value    :: numptr
    end function amg_zumf_solve
  end interface

  interface 
    function amg_zumf_free(symptr, numptr)&
         & bind(c,name='amg_zumf_free') result(info)
      use iso_c_binding
      integer(c_int)        :: info
      type(c_ptr), value    :: symptr, numptr
    end function amg_zumf_free
  end interface

  interface
    subroutine amg_z_umf_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      use psb_base_mod
      import amg_z_umf_solver_type
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_z_umf_solver_type), intent(inout) :: sv
      complex(psb_dpk_),intent(inout)         :: x(:)
      complex(psb_dpk_),intent(inout)         :: y(:)
      complex(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      complex(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      complex(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_z_umf_solver_apply
  end interface
  
  interface
    subroutine amg_z_umf_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      use psb_base_mod
      import amg_z_umf_solver_type
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_z_umf_solver_type), intent(inout) :: sv
      type(psb_z_vect_type),intent(inout)  :: x
      type(psb_z_vect_type),intent(inout)  :: y
      complex(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      complex(psb_dpk_),target, intent(inout) :: work(:)
      type(psb_z_vect_type),intent(inout) :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_z_vect_type),intent(inout), optional   :: initu
    end subroutine amg_z_umf_solver_apply_vect
  end interface

  interface 
    subroutine amg_z_umf_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      use psb_base_mod
      import amg_z_umf_solver_type
      Implicit None
      
      ! Arguments
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(amg_z_umf_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                        :: info
      type(psb_zspmat_type), intent(in), target, optional :: b
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_z_umf_solver_bld
  end interface

contains


  subroutine z_umf_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(amg_z_umf_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='z_umf_solver_free'

    call psb_erractionsave(err_act)

    call sv%clear_data(info) 
    
    if (info /= psb_success_) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_umf_solver_free


  subroutine z_umf_solver_clear_data(sv,info)

    Implicit None

    ! Arguments
    class(amg_z_umf_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='z_umf_solver_clear_data'

    call psb_erractionsave(err_act)
    info = 0 
    if (c_associated(sv%symbolic).and.c_associated(sv%numeric)) then 
      info = amg_zumf_free(sv%symbolic,sv%numeric)
      
      if (info /= psb_success_) goto 9999
      sv%symbolic = c_null_ptr
      sv%numeric  = c_null_ptr
      sv%symbsize = 0
      sv%numsize  = 0
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_umf_solver_clear_data

  subroutine z_umf_solver_finalize(sv)

    Implicit None

    ! Arguments
    type(amg_z_umf_solver_type), intent(inout) :: sv
    integer :: info
    Integer :: err_act
    character(len=20)  :: name='z_umf_solver_finalize'

    call sv%free(info) 

    return
  
  end subroutine z_umf_solver_finalize

  subroutine z_umf_solver_descr(sv,info,iout,coarse,prefix)

    Implicit None

    ! Arguments
    class(amg_z_umf_solver_type), intent(in) :: sv
    integer, intent(out)                       :: info
    integer, intent(in), optional              :: iout
    logical, intent(in), optional              :: coarse
      character(len=*), intent(in), optional   :: prefix

    ! Local variables
    integer      :: err_act
    character(len=20), parameter :: name='amg_z_umf_solver_descr'
    integer :: iout_
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
    
    write(iout_,*) trim(prefix_), '  UMFPACK Sparse Factorization Solver. '

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine z_umf_solver_descr

  function z_umf_solver_sizeof(sv) result(val)

    implicit none 
    ! Arguments
    class(amg_z_umf_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer             :: i

    val = 2*psb_sizeof_lp 
    val = val + sv%symbsize
    val = val + sv%numsize
    return
  end function z_umf_solver_sizeof

  function z_umf_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "UMFPACK solver"
  end function z_umf_solver_get_fmt

  function z_umf_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_umf_
  end function z_umf_solver_get_id
#endif
end module amg_z_umf_solver

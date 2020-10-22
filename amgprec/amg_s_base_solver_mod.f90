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
! File: amg_s_base_solver_mod.f90
!
! Module: amg_s_base_solver_mod
!
!  This module defines: 
!  - the amg_s_base_solver_type data structure containing the
!    basic solver type acting on a subdomain
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the solver is correctly defined;
!  - printing a	description of the solver;
!  - deallocating the data structure.  
!

module amg_s_base_solver_mod

  use amg_base_prec_type
  use psb_base_mod, only : psb_sspmat_type, &
       & psb_s_vect_type, psb_s_base_vect_type, psb_s_base_sparse_mat, &
       & psb_spk_, psb_i_base_vect_type, psb_erractionsave, psb_error_handler
  !
  ! 
  ! Type: amg_T_base_solver_type.
  ! 
  !  It holds the local solver; it has no mandatory components. 
  !
  !  type  amg_T_base_solver_type
  !  end type amg_T_base_solver_type
  !
  !    build      -   Compute the actual contents of the smoother; includes
  !                   invocation of the build method on the solver component. 
  !    free       -   Release memory
  !    apply      -   Apply the smoother to a vector (or to an array); includes
  !                   invocation of the apply method on the solver component. 
  !    descr      -   Prints a description of the object.
  !    default    -   Set default values
  !    dump       -   Dump to file object contents
  !    set        -   Sets various parameters; when a request is unknown
  !                   it is passed to the smoother object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !    stringval  -   convert string to val for internal parms
  !    get_fmt    -   short string descriptor
  !    get_id     -   numeric id descriptro
  !    get_wrksz  -   How many workspace vector does apply_vect need
  !
  !

  type amg_s_base_solver_type
  contains
    procedure, pass(sv) :: apply_v => amg_s_base_solver_apply_vect
    procedure, pass(sv) :: apply_a => amg_s_base_solver_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sv) :: check => amg_s_base_solver_check
    procedure, pass(sv) :: dump  => amg_s_base_solver_dmp
    procedure, pass(sv) :: clone => amg_s_base_solver_clone
    procedure, pass(sv) :: clone_settings => amg_s_base_solver_clone_settings
    procedure, pass(sv) :: clear_data     => amg_s_base_solver_clear_data
    procedure, pass(sv) :: build => amg_s_base_solver_bld
    procedure, pass(sv) :: cnv   => amg_s_base_solver_cnv
    procedure, pass(sv) :: free  => amg_s_base_solver_free
    procedure, pass(sv) :: cseti => amg_s_base_solver_cseti
    procedure, pass(sv) :: csetc => amg_s_base_solver_csetc
    procedure, pass(sv) :: csetr => amg_s_base_solver_csetr
    generic, public     :: set   => cseti, csetc, csetr
    procedure, pass(sv) :: default => s_base_solver_default
    procedure, pass(sv) :: descr   => amg_s_base_solver_descr
    procedure, pass(sv) :: sizeof  => s_base_solver_sizeof
    procedure, pass(sv) :: get_nzeros => s_base_solver_get_nzeros
    procedure, nopass   :: get_wrksz => s_base_solver_get_wrksize
    procedure, nopass   :: stringval => amg_stringval
    procedure, nopass   :: get_fmt   => s_base_solver_get_fmt
    procedure, nopass   :: get_id    => s_base_solver_get_id
    procedure, nopass   :: is_iterative => s_base_solver_is_iterative
    procedure, pass(sv) :: is_global => s_base_solver_is_global
  end type amg_s_base_solver_type

  private :: s_base_solver_sizeof, s_base_solver_default,&
       &  s_base_solver_get_nzeros, s_base_solver_get_fmt, &
       &  s_base_solver_is_iterative, s_base_solver_get_id, &
       &  s_base_solver_get_wrksize, s_base_solver_is_global


  interface  
    subroutine amg_s_base_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
       & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
       & amg_s_base_solver_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)           :: desc_data
      class(amg_s_base_solver_type), intent(inout) :: sv
      real(psb_spk_),intent(inout)              :: x(:)
      real(psb_spk_),intent(inout)              :: y(:)
      real(psb_spk_),intent(in)                 :: alpha,beta
      character(len=1),intent(in)                :: trans
      real(psb_spk_),target, intent(inout)      :: work(:)
      integer(psb_ipk_), intent(out)             :: info
      character, intent(in), optional       :: init
      real(psb_spk_),intent(inout), optional :: initu(:)
    end subroutine amg_s_base_solver_apply
  end interface 
  
      
  interface 
    subroutine amg_s_base_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)              :: desc_data
      class(amg_s_base_solver_type), intent(inout) :: sv
      type(psb_s_vect_type),intent(inout)          :: x
      type(psb_s_vect_type),intent(inout)          :: y
      real(psb_spk_),intent(in)                     :: alpha,beta
      character(len=1),intent(in)                    :: trans
      real(psb_spk_),target, intent(inout)          :: work(:)
      type(psb_s_vect_type),intent(inout)            :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_s_vect_type),intent(inout), optional   :: initu
    end subroutine amg_s_base_solver_apply_vect
  end interface
  
  interface 
    subroutine amg_s_base_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
       & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
       & amg_s_base_solver_type, psb_ipk_, psb_i_base_vect_type      
      Implicit None
      
      ! Arguments
      type(psb_sspmat_type), intent(in), target             :: a
      Type(psb_desc_type), Intent(inout)                    :: desc_a 
      class(amg_s_base_solver_type), intent(inout)          :: sv
      integer(psb_ipk_), intent(out)                        :: info
      type(psb_sspmat_type), intent(in), target, optional   :: b
      class(psb_s_base_sparse_mat), intent(in), optional    :: amold
      class(psb_s_base_vect_type), intent(in), optional     :: vmold
      class(psb_i_base_vect_type), intent(in), optional     :: imold
    end subroutine amg_s_base_solver_bld
  end interface
  
  interface 
    subroutine amg_s_base_solver_cnv(sv,info,amold,vmold,imold)
      import :: psb_s_base_sparse_mat, psb_s_base_vect_type, psb_spk_, &
       & amg_s_base_solver_type, psb_ipk_, psb_i_base_vect_type      
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout)          :: sv
      integer(psb_ipk_), intent(out)                        :: info
      class(psb_s_base_sparse_mat), intent(in), optional    :: amold
      class(psb_s_base_vect_type), intent(in), optional     :: vmold
      class(psb_i_base_vect_type), intent(in), optional     :: imold
    end subroutine amg_s_base_solver_cnv
  end interface
  
  interface 
    subroutine amg_s_base_solver_check(sv,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout)   :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_s_base_solver_check
  end interface
  
  interface 
    subroutine amg_s_base_solver_cseti(sv,what,val,info,idx)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                   :: what 
      integer(psb_ipk_), intent(in)                  :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine amg_s_base_solver_cseti
  end interface
  
  interface 
    subroutine amg_s_base_solver_csetc(sv,what,val,info,idx)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, & 
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                   :: what 
      character(len=*), intent(in)                   :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine amg_s_base_solver_csetc
  end interface 
  
  interface 
    subroutine amg_s_base_solver_csetr(sv,what,val,info,idx)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_            
      Implicit None      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                   :: what 
      real(psb_spk_), intent(in)                      :: val
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), intent(in), optional        :: idx
    end subroutine amg_s_base_solver_csetr
  end interface
  
  interface
    subroutine amg_s_base_solver_free(sv,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_s_base_solver_free
  end interface
  
  interface
    subroutine amg_s_base_solver_descr(sv,info,iout,coarse)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), intent(in), optional     :: iout
      logical, intent(in), optional               :: coarse

    end subroutine amg_s_base_solver_descr
  end interface 
  
  interface 
    subroutine amg_s_base_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_      
      implicit none 
      class(amg_s_base_solver_type), intent(in) :: sv
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(in)               :: level
      integer(psb_ipk_), intent(out)              :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: solver, global_num
    end subroutine amg_s_base_solver_dmp
  end interface
 
  interface
    subroutine amg_s_base_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout)              :: sv
      class(amg_s_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_s_base_solver_clone
  end interface

  interface
    subroutine amg_s_base_solver_clone_settings(sv,svout,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout) :: sv
      class(amg_s_base_solver_type), intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_s_base_solver_clone_settings
  end interface

  interface
    subroutine amg_s_base_solver_clear_data(sv,info)
      import :: psb_desc_type, psb_sspmat_type,  psb_s_base_sparse_mat, &
           & psb_s_vect_type, psb_s_base_vect_type, psb_spk_, &
           & amg_s_base_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_s_base_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_s_base_solver_clear_data
  end interface

contains
  !
  ! Function returning the size of the data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function s_base_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(amg_s_base_solver_type), intent(in) :: sv
    integer(psb_epk_)                  :: val
    integer(psb_ipk_)             :: i
    val = 0

    return
  end function s_base_solver_sizeof

  function s_base_solver_get_nzeros(sv) result(val)
    implicit none 
    class(amg_s_base_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i
    val = 0
  end function s_base_solver_get_nzeros

  subroutine s_base_solver_default(sv) 
    implicit none 
    ! Arguments
    class(amg_s_base_solver_type), intent(inout) :: sv
    ! Do nothing for base version

    return
  end subroutine s_base_solver_default

  function s_base_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Base solver"
  end function s_base_solver_get_fmt

  !
  ! If this is true, then the solver needs a starting
  ! guess. Currently only handled in JAC smoother. 
  ! 
  function s_base_solver_is_iterative() result(val)
    implicit none 
    logical  :: val

    val = .false.
  end function s_base_solver_is_iterative
  !
  ! Is the solver acting globally? In most cases 
  ! not, SuperLU_Dist does, MUMPS can do either.
  ! 
  function s_base_solver_is_global(sv) result(val)
    implicit none 
    class(amg_s_base_solver_type), intent(in) :: sv
    logical  :: val

    val = .false.
  end function s_base_solver_is_global

  function s_base_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_f_none_
  end function s_base_solver_get_id

  function s_base_solver_get_wrksize() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = 0
  end function s_base_solver_get_wrksize

end module amg_s_base_solver_mod

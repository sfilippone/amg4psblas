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
! File: amg_z_base_smoother_mod.f90
!
! Module: amg_z_base_smoother_mod
!
!  This module defines: 
!  - the amg_z_base_smoother_type data structure containing the
!    smoother  and related data structures;
!
!  It contains routines for
!  - Building and applying; 
!  - checking if the smoother is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.
!
!  What is the difference between a smoother and a solver?
!  In the mathematics literature the two concepts are treated
!  essentially as synonymous, but here we are using them in a more
!  computer-science oriented fashion. In particular, a SMOOTHER object
!  contains a SOLVER object: the SOLVER operates locally within the
!  current process, whereas the SMOOTHER object accounts for (possible)
!  interactions between processes.
!  Some solvers (MUMPS and SuperLU_DIST) can also operate on the entire
!  distributed matrix, in which case the smoother object essentially
!  becomes transparent. 
!
module amg_z_base_smoother_mod

  use amg_z_base_solver_mod
  use psb_base_mod, only : psb_desc_type, psb_zspmat_type, psb_epk_,&
       & psb_z_vect_type, psb_z_base_vect_type, psb_z_base_sparse_mat, &
       & psb_dpk_, psb_i_base_vect_type, psb_erractionsave, psb_error_handler
  
  !
  !
  ! 
  ! Type: amg_T_base_smoother_type.
  ! 
  !  It holds the smoother a single level. Its only mandatory component is a solver
  !  object which holds a local solver; this decoupling allows to have the same solver
  !  e.g ILU to work with Jacobi with multiple sweeps as well as with any AS variant.
  !
  !  type  amg_T_base_smoother_type
  !    class(amg_T_base_solver_type), allocatable :: sv
  !  end type amg_T_base_smoother_type
  !
  !   Methods:  
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
  !                   it is passed to the solver object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !    stringval  -   convert string to val for internal parms
  !    get_fmt    -   short string descriptor
  !    get_id     -   numeric id descriptro
  !    get_wrksz  -   How many workspace vector does apply_vect need
  !
  !
  ! 

  type  amg_z_base_smoother_type
    class(amg_z_base_solver_type), allocatable :: sv
  contains
    procedure, pass(sm) :: apply_v => amg_z_base_smoother_apply_vect
    procedure, pass(sm) :: apply_a => amg_z_base_smoother_apply
    generic, public     :: apply => apply_a, apply_v
    procedure, pass(sm) :: check => amg_z_base_smoother_check
    procedure, pass(sm) :: dump  => amg_z_base_smoother_dmp
    procedure, pass(sm) :: clone => amg_z_base_smoother_clone
    procedure, pass(sm) :: build => amg_z_base_smoother_bld
    procedure, pass(sm) :: cnv   => amg_z_base_smoother_cnv
    procedure, pass(sm) :: free  => amg_z_base_smoother_free
    procedure, pass(sm) :: clone_settings => amg_z_base_smoother_clone_settings
    procedure, pass(sm) :: clear_data     => amg_z_base_smoother_clear_data
    procedure, pass(sm) :: cseti => amg_z_base_smoother_cseti
    procedure, pass(sm) :: csetc => amg_z_base_smoother_csetc
    procedure, pass(sm) :: csetr => amg_z_base_smoother_csetr
    generic, public     :: set   => cseti, csetc, csetr
    procedure, pass(sm) :: default => z_base_smoother_default
    procedure, pass(sm) :: descr =>   amg_z_base_smoother_descr
    procedure, pass(sm) :: sizeof =>  z_base_smoother_sizeof
    procedure, pass(sm) :: get_nzeros => z_base_smoother_get_nzeros
    procedure, pass(sm) :: get_wrksz => z_base_smoother_get_wrksize
    procedure, nopass   :: stringval => amg_stringval
    procedure, nopass   :: get_fmt   => z_base_smoother_get_fmt
    procedure, nopass   :: get_id    => z_base_smoother_get_id
  end type amg_z_base_smoother_type


  private :: z_base_smoother_sizeof, z_base_smoother_get_fmt, &
       &  z_base_smoother_default, z_base_smoother_get_nzeros, &
       & z_base_smoother_get_id, z_base_smoother_get_wrksize



  interface 
    subroutine amg_z_base_smoother_apply(alpha,sm,x,beta,y,desc_data,& 
         & trans,sweeps,work,info,init,initu)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      type(psb_desc_type), intent(in)             :: desc_data
      class(amg_z_base_smoother_type), intent(inout) :: sm
      complex(psb_dpk_),intent(inout)                :: x(:)
      complex(psb_dpk_),intent(inout)                :: y(:)
      complex(psb_dpk_),intent(in)                   :: alpha,beta
      character(len=1),intent(in)                  :: trans
      integer(psb_ipk_), intent(in)                :: sweeps
      complex(psb_dpk_),target, intent(inout)        :: work(:)
      integer(psb_ipk_), intent(out)               :: info
      character, intent(in), optional       :: init
      complex(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_z_base_smoother_apply
  end interface
  
  interface 
    subroutine amg_z_base_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,&
         &  trans,sweeps,work,wv,info,init,initu)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      type(psb_desc_type), intent(in)                :: desc_data
      class(amg_z_base_smoother_type), intent(inout) :: sm
      type(psb_z_vect_type),intent(inout)            :: x
      type(psb_z_vect_type),intent(inout)            :: y
      complex(psb_dpk_),intent(in)                       :: alpha,beta
      character(len=1),intent(in)                      :: trans
      integer(psb_ipk_), intent(in)                    :: sweeps
      complex(psb_dpk_),target, intent(inout)            :: work(:)
      type(psb_z_vect_type),intent(inout)            :: wv(:)
      integer(psb_ipk_), intent(out)                   :: info
      character, intent(in), optional                :: init
      type(psb_z_vect_type),intent(inout), optional   :: initu
    end subroutine amg_z_base_smoother_apply_vect
  end interface
  
  interface 
    subroutine amg_z_base_smoother_check(sm,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine amg_z_base_smoother_check
  end interface
  
  interface 
    subroutine amg_z_base_smoother_cseti(sm,what,val,info,idx)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm 
      character(len=*), intent(in)                     :: what 
      integer(psb_ipk_), intent(in)                    :: val
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_), intent(in), optional          :: idx
    end subroutine amg_z_base_smoother_cseti
  end interface
  
  interface 
    subroutine amg_z_base_smoother_csetc(sm,what,val,info,idx)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      class(amg_z_base_smoother_type), intent(inout) :: sm 
      character(len=*), intent(in)                     :: what 
      character(len=*), intent(in)                     :: val
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_), intent(in), optional          :: idx
    end subroutine amg_z_base_smoother_csetc
  end interface
  
  interface 
    subroutine amg_z_base_smoother_csetr(sm,what,val,info,idx)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm 
      character(len=*), intent(in)                     :: what 
      real(psb_dpk_), intent(in)                        :: val
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_), intent(in), optional          :: idx
    end subroutine amg_z_base_smoother_csetr
  end interface
  
  interface 
    subroutine amg_z_base_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_, psb_i_base_vect_type
      ! Arguments
      type(psb_zspmat_type), intent(in), target     :: a
      Type(psb_desc_type), Intent(inout)              :: desc_a 
      class(amg_z_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(out)                  :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_z_base_smoother_bld
  end interface
  
  interface 
    subroutine amg_z_base_smoother_cnv(sm,info,amold,vmold,imold)
      import :: psb_z_base_sparse_mat, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_, psb_i_base_vect_type
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm 
      integer(psb_ipk_), intent(out)                  :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_z_base_smoother_cnv
  end interface
  
  interface 
    subroutine amg_z_base_smoother_free(sm,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine amg_z_base_smoother_free
  end interface
  
  interface 
    subroutine amg_z_base_smoother_descr(sm,info,iout,coarse)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      ! Arguments
      class(amg_z_base_smoother_type), intent(in) :: sm
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), intent(in), optional       :: iout
      logical, intent(in), optional                 :: coarse
    end subroutine amg_z_base_smoother_descr
  end interface
  
  interface 
    subroutine amg_z_base_smoother_dmp(sm,desc,level,info,prefix,head,smoother,solver,global_num)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      class(amg_z_base_smoother_type), intent(in) :: sm
      type(psb_desc_type), intent(in)               :: desc
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: smoother, solver, global_num
    end subroutine amg_z_base_smoother_dmp
  end interface
   
  interface
    subroutine amg_z_base_smoother_clone(sm,smout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout)              :: sm
      class(amg_z_base_smoother_type), allocatable, intent(inout) :: smout
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_z_base_smoother_clone
  end interface
   
  interface
    subroutine amg_z_base_smoother_clone_settings(sm,smout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm
      class(amg_z_base_smoother_type), intent(inout) :: smout
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine amg_z_base_smoother_clone_settings
  end interface
   
  interface
    subroutine amg_z_base_smoother_clear_data(sm,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, &
           & amg_z_base_smoother_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_z_base_smoother_type), intent(inout) :: sm
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine amg_z_base_smoother_clear_data
  end interface
  
contains
  !
  ! Function returning the size of the amg_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved. 
  !

  function z_base_smoother_get_nzeros(sm) result(val)
    implicit none 
    class(amg_z_base_smoother_type), intent(in) :: sm
    integer(psb_epk_) :: val
    integer(psb_ipk_)             :: i
    val = 0
    if (allocated(sm%sv)) &
         &  val =  sm%sv%get_nzeros()
  end function z_base_smoother_get_nzeros

  function z_base_smoother_sizeof(sm) result(val)
    implicit none 
    ! Arguments
    class(amg_z_base_smoother_type), intent(in) :: sm
    integer(psb_epk_)                    :: val
    integer(psb_ipk_)             :: i
    
    val = 0
    if (allocated(sm%sv)) then 
      val = sm%sv%sizeof()
    end if

    return
  end function z_base_smoother_sizeof

  !
  ! Set sensible defaults.
  ! To be called immediately after allocation
  !
  subroutine z_base_smoother_default(sm) 
    implicit none 
    ! Arguments
    class(amg_z_base_smoother_type), intent(inout) :: sm
    ! Do nothing for base version

    if (allocated(sm%sv)) call sm%sv%default()

    return
  end subroutine z_base_smoother_default

  function z_base_smoother_get_wrksize(sm) result(val)
    implicit none 
    class(amg_z_base_smoother_type), intent(inout) :: sm
    integer(psb_ipk_)  :: val

    val = 0
    if (allocated(sm%sv)) val = val + sm%sv%get_wrksz()
    
  end function z_base_smoother_get_wrksize
  
  function z_base_smoother_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Base smoother"
  end function z_base_smoother_get_fmt

  function z_base_smoother_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = amg_base_smooth_
  end function z_base_smoother_get_id

end module amg_z_base_smoother_mod

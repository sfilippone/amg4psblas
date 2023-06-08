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
!
!
module amg_c_ainv_solver

  use amg_c_base_ainv_mod
  use psb_base_mod, only : psb_d_vect_type

  type, extends(amg_c_base_ainv_solver_type) :: amg_c_ainv_solver_type
    !
    !  Compute an approximate factorization
    !      A^-1 = Z D^-1 W^T
    !  Note that here W is going to be transposed explicitly,
    !  so that the component w will in the end contain W^T.
    !
    integer(psb_ipk_)   :: alg, fill_in
    real(psb_spk_)      :: thresh
  contains
    procedure, pass(sv) :: check   => amg_c_ainv_solver_check
    procedure, pass(sv) :: build   => amg_c_ainv_solver_bld
    procedure, pass(sv) :: clone   => amg_c_ainv_solver_clone
    procedure, pass(sv) :: cseti   => amg_c_ainv_solver_cseti
    procedure, pass(sv) :: csetc   => amg_c_ainv_solver_csetc
    procedure, pass(sv) :: csetr   => amg_c_ainv_solver_csetr
!!$    procedure, pass(sv) :: seti    => amg_c_ainv_solver_seti
!!$    procedure, pass(sv) :: setc    => amg_c_ainv_solver_setc
!!$    procedure, pass(sv) :: setr    => amg_c_ainv_solver_setr
    procedure, pass(sv) :: descr   => amg_c_ainv_solver_descr
    procedure, pass(sv) :: default => c_ainv_solver_default
    procedure, nopass   :: stringval  => c_ainv_stringval
    procedure, nopass   :: algname => c_ainv_algname
  end type amg_c_ainv_solver_type


  private :: c_ainv_stringval, c_ainv_solver_default, &
       &  c_ainv_algname

  interface
    subroutine amg_c_ainv_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
       & amg_c_base_solver_type, psb_dpk_, amg_c_ainv_solver_type, psb_ipk_
      Implicit None
      class(amg_c_ainv_solver_type), intent(inout)              :: sv
      class(amg_c_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)                            :: info
    end subroutine amg_c_ainv_solver_clone
  end interface


  interface
    subroutine amg_c_ainv_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
       & psb_d_vect_type, psb_c_base_vect_type, psb_dpk_,&
       & amg_c_ainv_solver_type, psb_i_base_vect_type, psb_ipk_

      Implicit None

      ! Arguments
      type(psb_cspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a
      class(amg_c_ainv_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_cspmat_type), intent(in), target, optional :: b
      class(psb_c_base_sparse_mat), intent(in), optional  :: amold
      class(psb_c_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_c_ainv_solver_bld
  end interface

  interface
    subroutine amg_c_ainv_solver_check(sv,info)
      import :: psb_dpk_, amg_c_ainv_solver_type, psb_ipk_

      Implicit None

      ! Arguments
      class(amg_c_ainv_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_c_ainv_solver_check
  end interface

  interface
    subroutine amg_c_ainv_solver_cseti(sv,what,val,info,idx)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, psb_ipk_,&
           & psb_d_vect_type, psb_c_base_vect_type, psb_dpk_, amg_c_ainv_solver_type
      Implicit None
      ! Arguments
      class(amg_c_ainv_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                 :: what
      integer(psb_ipk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine amg_c_ainv_solver_cseti
  end interface


  interface
    subroutine amg_c_ainv_solver_csetc(sv,what,val,info,idx)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, psb_ipk_,&
           & psb_d_vect_type, psb_c_base_vect_type, psb_dpk_, amg_c_ainv_solver_type
      Implicit None
      ! Arguments
      class(amg_c_ainv_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                 :: what
      character(len=*), intent(in)                 :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine amg_c_ainv_solver_csetc
  end interface

  interface
    subroutine amg_c_ainv_solver_csetr(sv,what,val,info,idx)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat,  psb_ipk_,&
           & psb_d_vect_type, psb_c_base_vect_type, psb_spk_, amg_c_ainv_solver_type
      Implicit None
      ! Arguments
      class(amg_c_ainv_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                 :: what
      real(psb_spk_), intent(in)                   :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine amg_c_ainv_solver_csetr
  end interface

!!$  interface
!!$    subroutine amg_c_ainv_solver_setc(sv,what,val,info)
!!$      import :: amg_c_ainv_solver_type, psb_ipk_
!!$      Implicit none
!!$      ! Arguments
!!$      class(amg_c_ainv_solver_type), intent(inout) :: sv
!!$      integer(psb_ipk_), intent(in)                :: what
!!$      character(len=*), intent(in)                 :: val
!!$      integer(psb_ipk_), intent(out)               :: info
!!$    end subroutine amg_c_ainv_solver_setc
!!$  end interface
!!$
!!$  interface
!!$    subroutine amg_c_ainv_solver_seti(sv,what,val,info)
!!$      import :: amg_c_ainv_solver_type, psb_ipk_
!!$      Implicit none
!!$      ! Arguments
!!$      class(amg_c_ainv_solver_type), intent(inout) :: sv
!!$      integer(psb_ipk_), intent(in)                :: what
!!$      integer(psb_ipk_), intent(in)                :: val
!!$      integer(psb_ipk_), intent(out)               :: info
!!$    end subroutine amg_c_ainv_solver_seti
!!$  end interface
!!$
!!$  interface
!!$    subroutine amg_c_ainv_solver_setr(sv,what,val,info)
!!$      import :: amg_c_ainv_solver_type, psb_ipk_, psb_spk_
!!$      Implicit none
!!$      ! Arguments
!!$      class(amg_c_ainv_solver_type), intent(inout) :: sv
!!$      integer(psb_ipk_), intent(in)                :: what
!!$      real(psb_spk_), intent(in)                   :: val
!!$      integer(psb_ipk_), intent(out)               :: info
!!$    end subroutine amg_c_ainv_solver_setr
!!$  end interface

  interface
    subroutine amg_c_ainv_solver_descr(sv,info,iout,coarse,prefix)
      import :: psb_dpk_, amg_c_ainv_solver_type, psb_ipk_

      Implicit None

      ! Arguments
      class(amg_c_ainv_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: iout
      logical, intent(in), optional             :: coarse
      character(len=*), intent(in), optional  :: prefix
    end subroutine amg_c_ainv_solver_descr
  end interface

  interface  amg_ainv_bld
    subroutine amg_c_ainv_bld(a,alg,fillin,thresh,wmat,d,zmat,desc,info,blck,iscale)
      import :: psb_desc_type, psb_cspmat_type,  psb_c_base_sparse_mat, &
           & psb_d_vect_type, psb_c_base_vect_type, psb_spk_, psb_ipk_
      implicit none
      type(psb_cspmat_type), intent(in), target   :: a
      integer(psb_ipk_), intent(in)               :: fillin,alg
      real(psb_spk_), intent(in)                  :: thresh
      type(psb_cspmat_type), intent(inout)        :: wmat, zmat
      complex(psb_spk_), allocatable                 :: d(:)
      Type(psb_desc_type), Intent(inout)          :: desc
      integer(psb_ipk_), intent(out)              :: info
      type(psb_cspmat_type), intent(in), optional :: blck
      integer(psb_ipk_), intent(in), optional     :: iscale
    end subroutine amg_c_ainv_bld
  end interface


contains

  subroutine c_ainv_solver_default(sv)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(amg_c_ainv_solver_type), intent(inout) :: sv

    sv%alg     = amg_ainv_llk_
    sv%fill_in = 0
    sv%thresh  = dzero

    return
  end subroutine c_ainv_solver_default

  function is_positive_nz_min(ip) result(res)
    implicit none
    integer(psb_ipk_), intent(in) :: ip
    logical             :: res

    res = (ip >= 1)
    return
  end function is_positive_nz_min


  function c_ainv_stringval(string) result(val)
    use psb_base_mod, only : psb_ipk_,psb_toupper
    implicit none
  ! Arguments
    character(len=*), intent(in) :: string
    integer(psb_ipk_) :: val
    character(len=*), parameter :: name='d_ainv_stringval'

    select case(psb_toupper(trim(string)))
    case('LLK')
      val = amg_ainv_llk_
    case('STAB-LLK')
      val = amg_ainv_s_ft_llk_
    case('SYM-LLK')
      val = amg_ainv_s_llk_
    case('MLK')
      val = amg_ainv_mlk_
#if defined(HAVE_TUMA_SAINV)
    case('SAINV-TUMA')
      val = amg_ainv_s_tuma_
    case('LAINV-TUMA')
      val = amg_ainv_l_tuma_
#endif
    case default
      val  = amg_stringval(string)
    end select
  end function c_ainv_stringval


  function c_ainv_algname(ialg) result(val)
    integer(psb_ipk_), intent(in) :: ialg
    character(len=40) :: val

    character(len=*), parameter :: mlkname   = 'Left-looking, list merge '
    character(len=*), parameter :: llkname   = 'Left-looking '
    character(len=*), parameter :: stabllkname  = 'Stabilized Left-looking '
    character(len=*), parameter :: sllkname  = 'Symmetric Left-looking '
    character(len=*), parameter :: sainvname = 'SAINV (Benzi & Tuma) '
    character(len=*), parameter :: lainvname = 'LAINV (Benzi & Tuma) '
    character(len=*), parameter :: defname   = 'Unknown alg variant '

    select case (ialg)
    case(amg_ainv_mlk_)
      val = mlkname
    case(amg_ainv_llk_)
      val = llkname
    case(amg_ainv_s_llk_)
      val = sllkname
    case(amg_ainv_s_ft_llk_)
      val = stabllkname
#if defined(HAVE_TUMA_SAINV)
    case(amg_ainv_s_tuma_ )
      val = sainvname
    case(amg_ainv_l_tuma_ )
      val = lainvname
#endif
    case default
      val = defname
    end select

  end function c_ainv_algname

end module amg_c_ainv_solver

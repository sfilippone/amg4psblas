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

module amg_z_invt_solver

  use amg_z_base_solver_mod
  use amg_z_base_ainv_mod
  use psb_base_mod, only : psb_z_vect_type

  type, extends(amg_z_base_ainv_solver_type) :: amg_z_invt_solver_type
    integer(psb_ipk_)   :: fill_in, inv_fill
    real(psb_dpk_)      :: thresh, inv_thresh
  contains
    procedure, pass(sv) :: check   => amg_z_invt_solver_check
    procedure, pass(sv) :: clone   => amg_z_invt_solver_clone
    procedure, pass(sv) :: build   => amg_z_invt_solver_bld
    procedure, pass(sv) :: cseti   => amg_z_invt_solver_cseti
    procedure, pass(sv) :: csetr   => amg_z_invt_solver_csetr
    procedure, pass(sv) :: descr   => amg_z_invt_solver_descr
    procedure, pass(sv) :: default => z_invt_solver_default
  end type amg_z_invt_solver_type

  private :: z_invt_solver_default


  interface
    subroutine amg_z_invt_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & amg_z_base_solver_type, psb_dpk_, amg_z_invt_solver_type, psb_ipk_
      Implicit None
      class(amg_z_invt_solver_type), intent(inout)              :: sv
      class(amg_z_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)                            :: info
    end subroutine amg_z_invt_solver_clone
  end interface

  interface
    subroutine amg_z_invt_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat, &
       & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_,&
       & amg_z_invt_solver_type, psb_i_base_vect_type, psb_ipk_

      Implicit None

      ! Arguments
      type(psb_zspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a
      class(amg_z_invt_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_zspmat_type), intent(in), target, optional :: b
      class(psb_z_base_sparse_mat), intent(in), optional  :: amold
      class(psb_z_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_z_invt_solver_bld
  end interface

  interface
    subroutine amg_z_invt_solver_check(sv,info)
      import :: amg_z_invt_solver_type, psb_ipk_

      Implicit None

      ! Arguments
      class(amg_z_invt_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_z_invt_solver_check
  end interface

  interface
    subroutine amg_z_invt_solver_cseti(sv,what,val,info,idx)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat,  psb_ipk_,&
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, amg_z_invt_solver_type
      Implicit None
      ! Arguments
      class(amg_z_invt_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                 :: what
      integer(psb_ipk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine amg_z_invt_solver_cseti
  end interface

  interface
    subroutine amg_z_invt_solver_csetr(sv,what,val,info,idx)
      import :: psb_desc_type, psb_zspmat_type,  psb_z_base_sparse_mat,  psb_ipk_,&
           & psb_z_vect_type, psb_z_base_vect_type, psb_dpk_, amg_z_invt_solver_type
      Implicit None
      ! Arguments
      class(amg_z_invt_solver_type), intent(inout) :: sv
      character(len=*), intent(in)                 :: what
      real(psb_dpk_), intent(in)                   :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine amg_z_invt_solver_csetr
  end interface

  interface
    subroutine amg_z_invt_solver_descr(sv,info,iout,coarse,prefix)
      import :: psb_dpk_, amg_z_invt_solver_type, psb_ipk_

      Implicit None

      ! Arguments
      class(amg_z_invt_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), intent(in), optional     :: iout
      logical, intent(in), optional               :: coarse
      character(len=*), intent(in), optional      :: prefix

    end subroutine amg_z_invt_solver_descr
  end interface

contains

  subroutine z_invt_solver_default(sv)

    !use psb_base_mod

    Implicit None

    ! Arguments
    class(amg_z_invt_solver_type), intent(inout) :: sv

    sv%fill_in    = 0
    sv%inv_fill   = 0
    sv%thresh     = dzero
    sv%inv_thresh = dzero

    return
  end subroutine z_invt_solver_default

end module amg_z_invt_solver

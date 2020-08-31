!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
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
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
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
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. This variant is based on the hybrid method
!  presented in 
!
!    S. Gratton, P. Henon, P. Jiranek and X. Vasseur:
!    Reducing complexity of algebraic multigrid by aggregation
!    Numerical Lin. Algebra with Applications, 2016, 23:501-518
!    
module amg_z_hybrid_aggregator_mod

  use amg_z_dec_aggregator_mod
  !
  !   sm           -  class(amg_T_base_smoother_type), allocatable
  !                   The current level preconditioner (aka smoother).
  !   parms        -  type(amg_RTml_parms)
  !                   The parameters defining the multilevel strategy.
  !   ac           -  The local part of the current-level matrix, built by
  !                   coarsening the previous-level matrix.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the matrix
  !                   stored in ac.
  !   base_a       -  type(psb_Tspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the current 
  !                   matrix (so we have a unified treatment of residuals).
  !                   We need this to avoid passing explicitly the current matrix
  !                   to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the
  !                   matrix pointed by base_a.
  !   map          -  Stores the maps (restriction and prolongation) between the
  !                   vector spaces associated to the index spaces of the previous
  !                   and current levels.
  !
  !   Methods:  
  !     Most methods follow the encapsulation hierarchy: they take whatever action
  !     is appropriate for the current object, then call the corresponding method for
  !     the contained object.
  !     As an example: the descr() method prints out a description of the
  !     level. It starts by invoking the descr() method of the parms object,
  !     then calls the descr() method of the smoother object. 
  !
  !    descr      -   Prints a description of the object.
  !    default    -   Set default values
  !    dump       -   Dump to file object contents
  !    set        -   Sets various parameters; when a request is unknown
  !                   it is passed to the smoother object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !
  !
  type, extends(amg_z_dec_aggregator_type) :: amg_z_hybrid_aggregator_type
    
  contains
    procedure, pass(ag) :: bld_tprol    => amg_z_hybrid_aggregator_build_tprol
    procedure, nopass   :: fmt          => amg_z_hybrid_aggregator_fmt
  end type amg_z_hybrid_aggregator_type


  interface
    subroutine  amg_z_hybrid_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: amg_z_hybrid_aggregator_type, psb_desc_type, psb_zspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, amg_dml_parms
      implicit none
      class(amg_z_hybrid_aggregator_type), target, intent(inout) :: ag
      type(amg_dml_parms), intent(inout)  :: parms 
      type(psb_zspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_zspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine amg_z_hybrid_aggregator_build_tprol
  end interface

contains


  function amg_z_hybrid_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Hybrid Decoupled aggregation"
  end function amg_z_hybrid_aggregator_fmt

  
end module amg_z_hybrid_aggregator_mod

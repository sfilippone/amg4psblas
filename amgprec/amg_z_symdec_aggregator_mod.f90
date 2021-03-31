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
!  Locally symmetrized (decoupled) aggregation algorithm.
!  This version differs from the basic decoupled aggregation algorithm
!  only because it works on (the pattern of)  A+A^T  instead of A. 
!
!    
module amg_z_symdec_aggregator_mod

  use amg_z_dec_aggregator_mod
  !> \namespace  amg_z_symdec_aggregator_mod  \class  amg_z_symdec_aggregator_type
  !! \extends amg_z_dec_aggregator_mod::amg_z_dec_aggregator_type
  !!
  !!  This version differs from the basic decoupled aggregation algorithm
  !!   only because it works on (the pattern of)  A+A^T  instead of A.
  !!   
  !
  type, extends(amg_z_dec_aggregator_type) :: amg_z_symdec_aggregator_type
    
  contains
    procedure, pass(ag) :: bld_tprol    => amg_z_symdec_aggregator_build_tprol
    procedure, pass(ag) :: descr        => amg_z_symdec_aggregator_descr
    procedure, nopass   :: fmt          => amg_z_symdec_aggregator_fmt
  end type amg_z_symdec_aggregator_type


  interface
    subroutine  amg_z_symdec_aggregator_build_tprol(ag,parms,ag_data,&
         & a,desc_a,ilaggr,nlaggr,t_prol,info)
      import :: amg_z_symdec_aggregator_type, psb_desc_type, psb_zspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_lpk_,  psb_lzspmat_type, amg_dml_parms, amg_daggr_data
      implicit none
      class(amg_z_symdec_aggregator_type), target, intent(inout) :: ag
      type(amg_dml_parms), intent(inout)  :: parms 
      type(amg_daggr_data), intent(in)    :: ag_data
      type(psb_zspmat_type), intent(inout) :: a
      type(psb_desc_type), intent(inout)     :: desc_a
      integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_lzspmat_type), intent(out)  :: t_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine amg_z_symdec_aggregator_build_tprol
  end interface


contains

  function amg_z_symdec_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Symmetric Decoupled aggregation"
  end function amg_z_symdec_aggregator_fmt
  
  subroutine  amg_z_symdec_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(amg_z_symdec_aggregator_type), intent(in) :: ag
    type(amg_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'Decoupled Aggregator locally-symmetrized'
    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine amg_z_symdec_aggregator_descr

end module amg_z_symdec_aggregator_mod

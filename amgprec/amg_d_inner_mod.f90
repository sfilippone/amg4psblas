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
! File: amg_inner_mod.f90
!
! Module: amg_inner_mod
!
!  This module defines the interfaces to inner AMG4PSBLAS routines.
!  The interfaces  of the user level routines are defined in amg_prec_mod.f90.
!
module amg_d_inner_mod

  use psb_base_mod, only : psb_dspmat_type, psb_desc_type, psb_i_base_vect_type, &
       & psb_dpk_, psb_d_base_sparse_mat, psb_d_base_vect_type, psb_ipk_, &
       & psb_d_vect_type, psb_lpk_, psb_ldspmat_type
  use amg_d_prec_type, only : amg_dprec_type, amg_dml_parms, &
       & amg_d_onelev_type, amg_dmlprec_wrk_type

  interface amg_mlprec_bld
    subroutine amg_dmlprec_bld(a,desc_a,prec,info, amold, vmold,imold)
      import :: psb_dspmat_type, psb_desc_type, psb_i_base_vect_type, &
           & psb_dpk_, psb_d_base_sparse_mat, psb_d_base_vect_type, psb_ipk_
      import :: amg_dprec_type
      implicit none
      type(psb_dspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      type(amg_dprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine amg_dmlprec_bld
  end interface amg_mlprec_bld

  interface amg_mlprec_aply
    subroutine amg_dmlprec_aply(alpha,p,x,beta,y,desc_data,trans,work,info)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, psb_ipk_
      import :: amg_dprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(amg_dprec_type), intent(inout) :: p
      real(psb_dpk_),intent(in)         :: alpha,beta
      real(psb_dpk_),intent(inout)      :: x(:)
      real(psb_dpk_),intent(inout)      :: y(:)
      character,intent(in)               :: trans
      real(psb_dpk_),target             :: work(:)
      integer(psb_ipk_), intent(out)     :: info
    end subroutine amg_dmlprec_aply
    subroutine amg_dmlprec_aply_vect(alpha,p,x,beta,y,desc_data,trans,work,info)
      import :: psb_dspmat_type, psb_desc_type, &
           & psb_dpk_, psb_d_vect_type, psb_ipk_
      import :: amg_dprec_type
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      type(amg_dprec_type), intent(inout) :: p
      real(psb_dpk_),intent(in)            :: alpha,beta
      type(psb_d_vect_type),intent(inout) :: x
      type(psb_d_vect_type),intent(inout) :: y
      character,intent(in)                  :: trans
      real(psb_dpk_),target                :: work(:)
      integer(psb_ipk_), intent(out)        :: info
    end subroutine amg_dmlprec_aply_vect
  end interface amg_mlprec_aply
 
  interface amg_map_to_tprol
    subroutine amg_d_map_to_tprol(desc_a,ilaggr,nlaggr,op_prol,info)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, psb_ipk_, psb_lpk_, psb_ldspmat_type
      import :: amg_d_onelev_type
      implicit none 
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_lpk_), allocatable, intent(inout) :: ilaggr(:),nlaggr(:)
      type(psb_ldspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine amg_d_map_to_tprol
  end interface amg_map_to_tprol

  abstract interface
    subroutine amg_daggrmat_var_bld(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: psb_dspmat_type, psb_desc_type, psb_dpk_, psb_ipk_, psb_lpk_, psb_ldspmat_type
      import ::  amg_d_onelev_type, amg_dml_parms
      implicit none 
      type(psb_dspmat_type), intent(in)         :: a
      type(psb_desc_type), intent(inout)          :: desc_a
      integer(psb_lpk_), intent(inout)            :: ilaggr(:), nlaggr(:)
      type(amg_dml_parms), intent(inout)        :: parms 
      type(psb_dspmat_type), intent(inout)     :: op_prol,ac,op_restr
      type(psb_ldspmat_type), intent(inout)     :: t_prol
      type(psb_desc_type), intent(inout)         :: desc_ac
      integer(psb_ipk_), intent(out)             :: info
    end subroutine amg_daggrmat_var_bld
  end interface

  procedure(amg_daggrmat_var_bld) ::  amg_daggrmat_nosmth_bld, &
       & amg_daggrmat_smth_bld, amg_daggrmat_minnrg_bld

end module amg_d_inner_mod

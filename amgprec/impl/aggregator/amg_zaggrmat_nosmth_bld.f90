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
! File: amg_zaggrmat_nosmth_bld.F90
!
! Subroutine: amg_zaggrmat_nosmth_bld
! Version:    complex
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is the piecewise constant interpolation operator corresponding
!  the fine-to-coarse level mapping built by amg_aggrmap_bld.
! 
!  The coarse-level matrix A_C is distributed among the parallel processes or
!  replicated on each of them, according to the value of p%parms%coarse_mat
!  specified by the user through amg_zprecinit and amg_zprecset.
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  aggregator%mat_bld.
!
!  For details see
!    P. D'Ambra, D. di Serafino and  S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.,
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    a          -  type(psb_zspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(amg_z_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(amg_dml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_zspmat_type), output
!                  The coarse matrix on output 
!                  
!    ilaggr     -  integer, dimension(:), input
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix. Note that the indices
!                  are assumed to be shifted so as to make sure the ranges on
!                  the various processes do not   overlap.
!    nlaggr     -  integer, dimension(:) input
!                  nlaggr(i) contains the aggregates held by process i.
!    op_prol    -  type(psb_zspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_zspmat_type), output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!
!
subroutine amg_zaggrmat_nosmth_bld(a,desc_a,ilaggr,nlaggr,parms,&
     & ac,desc_ac,op_prol,op_restr,t_prol,info)
  use psb_base_mod
  use amg_base_prec_type
  use amg_z_inner_mod, amg_protect_name => amg_zaggrmat_nosmth_bld
  use amg_z_base_aggregator_mod
  implicit none

  ! Arguments
  type(psb_zspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(inout)         :: desc_a
  integer(psb_lpk_), intent(inout)           :: ilaggr(:), nlaggr(:)
  type(amg_dml_parms), intent(inout)      :: parms 
  type(psb_zspmat_type), intent(inout)       :: op_prol,ac,op_restr
  type(psb_lzspmat_type), intent(inout)    :: t_prol
  type(psb_desc_type), intent(inout)       :: desc_ac
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  integer(psb_ipk_)   :: err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  character(len=20)   :: name
  type(psb_lz_coo_sparse_mat) :: lcoo_prol
  type(psb_z_coo_sparse_mat) :: coo_prol, coo_restr
  type(psb_z_csr_sparse_mat) :: acsr
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nzl, ip, &
       & naggr, nzt, naggrm1, naggrp1, i, k
  integer(psb_ipk_) :: inaggr, nzlp
  logical, parameter :: debug = .false.

  name = 'amg_aggrmat_nosmth_bld'
  info = psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))

  call a%cp_to(acsr)
  call t_prol%mv_to(lcoo_prol)
  inaggr = naggr
  call psb_cdall(ctxt,desc_ac,info,nl=inaggr)
  nzlp = lcoo_prol%get_nzeros()
  call desc_ac%indxmap%g2lip_ins(lcoo_prol%ja(1:nzlp),info) 
  call lcoo_prol%set_ncols(desc_ac%get_local_cols())
  call lcoo_prol%cp_to_icoo(coo_prol,info)
  
  if (debug) call check_coo(me,trim(name)//' Check 1 on  coo_prol:',coo_prol)

  call psb_cdasb(desc_ac,info)
  call psb_cd_reinit(desc_ac,info)

  call amg_ptap_bld(acsr,desc_a,nlaggr,parms,ac,&
       & coo_prol,desc_ac,coo_restr,info)

  call coo_restr%set_nrows(desc_ac%get_local_rows())
  call coo_restr%set_ncols(desc_a%get_local_cols())
  call coo_prol%set_nrows(desc_a%get_local_rows())
  call coo_prol%set_ncols(desc_ac%get_local_cols())

  if (debug) call check_coo(me,trim(name)//' Check 1 on coo_restr:',coo_restr)

  call op_prol%mv_from(coo_prol)
  call op_restr%mv_from(coo_restr)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
contains
  subroutine check_coo(me,string,coo)
    implicit none
    integer(psb_ipk_) :: me
    type(psb_z_coo_sparse_mat) :: coo
    character(len=*) :: string
    integer(psb_lpk_) :: nr,nc,nz
    nr = coo%get_nrows()
    nc = coo%get_ncols()
    nz = coo%get_nzeros()
    write(0,*) me,string,nr,nc,&
         & minval(coo%ia(1:nz)),maxval(coo%ia(1:nz)),&
         & minval(coo%ja(1:nz)),maxval(coo%ja(1:nz))

  end subroutine check_coo
end subroutine amg_zaggrmat_nosmth_bld

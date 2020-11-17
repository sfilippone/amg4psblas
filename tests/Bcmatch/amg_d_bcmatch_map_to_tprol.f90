!   
!   
!                             MLD2P4  version 2.1
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 , 2017 
!  
!                        Salvatore Filippone  Cranfield University
!  		      Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
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
! File: mld_d_bcmatch_map_to_tprol.f90
!
! Subroutine: mld_d_bcmatch_map_to_tprol
! Version:    real
!
!  This routine uses a mapping from the row indices of the fine-level matrix
!  to the row indices of the coarse-level matrix to build a tentative
!  prolongator, i.e. a piecewise constant operator.
!  This is later used to build the final operator; the code has been refactored here
!  to  be shared among all the methods that provide the tentative prolongator
!  through a simple integer mapping.
!
!  The aggregation algorithm is a parallel version of that described in
!  * M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!  * P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
!    (1996), 179-196.
!  For more details see
!    P. D'Ambra, D. di Serafino and S. Filippone, On the development of
!    PSBLAS-based parallel two-level Schwarz preconditioners, Appl. Num. Math.
!    57 (2007), 1181-1196.
!
!
! Arguments:
!    aggr_type  -  integer, input.
!                  The scalar used to identify the aggregation algorithm.
!    theta      -  real, input.
!                  The aggregation threshold used in the aggregation algorithm.
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    ilaggr     -  integer, dimension(:), allocatable.
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix. Note that on exit the indices
!                  will be shifted so as to make sure the ranges on the various processes do not
!                  overlap.
!    nlaggr     -  integer, dimension(:), allocatable.
!                  nlaggr(i) contains the aggregates held by process i.
!    op_prol    -  type(psb_dspmat_type).
!               The tentative prolongator, based on ilaggr.
!               
!    info       -  integer, output.
!                  Error code.
!
subroutine mld_d_bcmatch_map_to_tprol(desc_a,ilaggr,nlaggr,valaggr, op_prol,info)

  use psb_base_mod
  use mld_d_inner_mod!, mld_protect_name => mld_d_bcmatch_map_to_tprol
  use mld_d_bcmatch_aggregator_mod, mld_protect_name => mld_d_bcmatch_map_to_tprol

  implicit none

  ! Arguments
  type(psb_desc_type), intent(in)    :: desc_a
  integer(psb_ipk_), allocatable, intent(inout)  :: ilaggr(:),nlaggr(:)
  real(psb_dpk_), allocatable, intent(inout)  :: valaggr(:)
  type(psb_dspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)               :: info

  ! Local variables
  integer(psb_ipk_) :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m,naggrm1, naggrp1, ntaggr
  type(psb_d_coo_sparse_mat) :: tmpcoo
  integer(psb_ipk_) :: debug_level, debug_unit,err_act
  integer(psb_ipk_) :: ctxt,np,me
  integer(psb_ipk_) :: nrow, ncol, n_ne
  character(len=20)  :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'mld_d_bcmatch_map_to_tprol'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ctxt=desc_a%get_context()
  call psb_info(ctxt,me,np)
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))
  ilaggr(1:nrow) = ilaggr(1:nrow) + naggrm1
  call psb_halo(ilaggr,desc_a,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
    goto 9999
  end if

  call psb_halo(valaggr,desc_a,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_halo')
    goto 9999
  end if

  call tmpcoo%allocate(ncol,ntaggr,ncol)
  do i=1,ncol
    tmpcoo%val(i) = valaggr(i)
    tmpcoo%ia(i)  = i
    tmpcoo%ja(i)  = ilaggr(i)  
  end do
  call tmpcoo%set_nzeros(ncol)
  call tmpcoo%set_dupl(psb_dupl_add_)
  call tmpcoo%set_sorted() ! At this point this is in row-major
  call op_prol%mv_from(tmpcoo)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine mld_d_bcmatch_map_to_tprol

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
! File: amg_z_map_to_tprol.f90
!
! Subroutine: amg_z_map_to_tprol
! Version:    complex
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
!    op_prol    -  type(psb_zspmat_type).
!               The tentative prolongator, based on ilaggr.
!               
!    info       -  integer, output.
!                  Error code.
!
subroutine amg_z_map_to_tprol(desc_a,ilaggr,nlaggr,op_prol,info)

  use psb_base_mod
  use amg_z_inner_mod, amg_protect_name => amg_z_map_to_tprol

  implicit none

  ! Arguments
  type(psb_desc_type), intent(in)    :: desc_a
  integer(psb_lpk_), allocatable, intent(inout)  :: ilaggr(:),nlaggr(:)
  type(psb_lzspmat_type), intent(out)  :: op_prol
  integer(psb_ipk_), intent(out)               :: info

  ! Local variables
  integer(psb_lpk_) :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m,naggrm1, naggrp1, ntaggr
  type(psb_lz_coo_sparse_mat) :: tmpcoo
  integer(psb_ipk_) :: debug_level, debug_unit,err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_lpk_)   :: nrow, ncol, n_ne
  character(len=20)   :: name, ch_err

  info=psb_success_
  name = 'amg_map_to_tprol'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
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

  call tmpcoo%allocate(nrow,ntaggr,ncol)
  k = 0 
  do i=1,nrow
    !
    ! Note: at this point, a value ilaggr(i)<=0
    ! tags a "singleton" row, and it has to be
    ! left alone. 
    !
    if (ilaggr(i)>0) then
      k = k + 1
      tmpcoo%val(k) = zone
      tmpcoo%ia(k)  = i
      tmpcoo%ja(k)  = ilaggr(i)
    end if
  end do
  call tmpcoo%set_nzeros(k)
  call tmpcoo%set_dupl(psb_dupl_add_)
  call tmpcoo%set_sorted() ! At this point this is in row-major
  call op_prol%mv_from(tmpcoo)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_z_map_to_tprol

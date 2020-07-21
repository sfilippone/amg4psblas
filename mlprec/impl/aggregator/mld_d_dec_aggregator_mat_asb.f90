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
! File: mld_d_dec_aggregator_mat_asb.f90
!
! Subroutine: mld_d_dec_aggregator_mat_asb
! Version:    real
!
!
!  From a given AC to final format, generating DESC_AC
! 
! Arguments:
!    ag       -  type(mld_d_dec_aggregator_type), input/output.
!               The aggregator object
!    parms   -  type(mld_dml_parms), input 
!               The aggregation parameters
!    a          -  type(psb_dspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
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
!    ac         -  type(psb_dspmat_type), inout
!                  The coarse matrix 
!    desc_ac    -  type(psb_desc_type), output.
!                  The communication descriptor of the fine-level matrix.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!                  
!    op_prol    -  type(psb_dspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_dspmat_type), input/output
!                  The restrictor operator; normally, it is the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!  
subroutine  mld_d_dec_aggregator_mat_asb(ag,parms,a,desc_a,&
     & ac,desc_ac, op_prol,op_restr,info)
  use psb_base_mod
  use mld_base_prec_type
  use mld_d_dec_aggregator_mod, mld_protect_name => mld_d_dec_aggregator_mat_asb
  implicit none
  class(mld_d_dec_aggregator_type), target, intent(inout) :: ag
  type(mld_dml_parms), intent(inout)  :: parms 
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_desc_type), intent(inout)     :: desc_a
  type(psb_dspmat_type), intent(inout) :: op_prol, ac,op_restr
  type(psb_desc_type), intent(inout)     :: desc_ac
  integer(psb_ipk_), intent(out)         :: info
  !
  integer(psb_ipk_)              :: ictxt, np, me
  type(psb_ld_coo_sparse_mat)  :: tmpcoo
  type(psb_ldspmat_type)       :: tmp_ac
  integer(psb_ipk_)              :: i_nr, i_nc, i_nl, nzl
  integer(psb_lpk_)              :: ntaggr
  integer(psb_ipk_) :: err_act, debug_level, debug_unit
  character(len=20) :: name='d_dec_aggregator_mat_asb'


  if (psb_get_errstatus().ne.0) return 
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  select case(parms%coarse_mat)

  case(mld_distr_mat_) 

    call ac%cscnv(info,type='csr')
    call op_prol%cscnv(info,type='csr')
    call op_restr%cscnv(info,type='csr')

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done ac '

  case(mld_repl_mat_) 
    !
    ! We are assuming here that an d matrix
    ! can hold all entries
    !
    if (desc_ac%get_global_rows() < huge(1_psb_ipk_) ) then 
      ntaggr = desc_ac%get_global_rows()
      i_nr   = ntaggr
    else
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
      goto 9999
    end if
    
    call op_prol%mv_to(tmpcoo)
    nzl = tmpcoo%get_nzeros()
    call psb_loc_to_glob(tmpcoo%ja(1:nzl),desc_ac,info,'I')
    call tmpcoo%set_ncols(i_nr)
    call op_prol%mv_from(tmpcoo)

    call op_restr%mv_to(tmpcoo)
    nzl = tmpcoo%get_nzeros()
    call psb_loc_to_glob(tmpcoo%ia(1:nzl),desc_ac,info,'I')
    call tmpcoo%set_nrows(i_nr)
    call op_restr%mv_from(tmpcoo)


    call psb_gather(tmp_ac,ac,desc_ac,info,root=-ione,&
         & dupl=psb_dupl_add_,keeploc=.false.)
    call tmp_ac%mv_to(tmpcoo)
    call ac%mv_from(tmpcoo)
    
    call psb_cdall(ictxt,desc_ac,info,mg=ntaggr,repl=.true.)
    if (info == psb_success_) call psb_cdasb(desc_ac,info)

    if (info /= psb_success_) goto 9999

  case default 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid mld_coarse_mat_')
    goto 9999
  end select


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return


contains
  subroutine check_coo(me,string,coo)
    implicit none
    integer(psb_ipk_) :: me
    type(psb_ld_coo_sparse_mat) :: coo
    character(len=*) :: string
    integer(psb_lpk_) :: nr,nc,nz
    nr = coo%get_nrows()
    nc = coo%get_ncols()
    nz = coo%get_nzeros()
    write(0,*) me,string,nr,nc,&
         & minval(coo%ia(1:nz)),maxval(coo%ia(1:nz)),&
         & minval(coo%ja(1:nz)),maxval(coo%ja(1:nz))

  end subroutine check_coo
end subroutine mld_d_dec_aggregator_mat_asb

!   
!   
!                             MLD2P4  Extensions
!    
!    (C) Copyright 2019
!  
!                        Salvatore Filippone  Cranfield University
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
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
! File: amg_rap.F90
!
!  This routine computes the triple product :
!                   AC = R  A P    
!   R is stored in Restr
!   WARNINGS:
!   1. So far this has only been tested with R=P^T
!   2. Patterns and descriptors need to be consistent
!      between initial build and application of this. 
!
subroutine amg_c_rap(a_csr,desc_a,nlaggr,parms,ac,&
     & coo_prol,desc_ac,coo_restr,info)
  use psb_base_mod
  use amg_c_inner_mod
  use amg_c_base_aggregator_mod, amg_protect_name => amg_c_rap
  implicit none

  ! Arguments
  type(psb_c_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(inout)          :: desc_a
  integer(psb_lpk_), intent(inout)           :: nlaggr(:)
  type(amg_sml_parms), intent(inout)         :: parms 
  type(psb_c_coo_sparse_mat), intent(inout)  :: coo_prol, coo_restr
  type(psb_desc_type), intent(inout)         :: desc_ac
  type(psb_cspmat_type), intent(out)        :: ac
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  integer(psb_ipk_)   :: err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, ndx
  character(len=40)   :: name
  type(psb_lc_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_c_csr_sparse_mat) :: acsr3, csr_prol, ac_csr, csr_restr
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nglob, ntaggr, naggrm1, naggrp1
  integer(psb_ipk_) :: nrow, ncol, nrl, nzl, ip, nzt, i, k
  integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza
  logical, parameter :: debug=.false.

  name='amg_rap'
  if(psb_get_errstatus().ne.0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1)) 
  !write(0,*)me,' ',name,' input sizes',nlaggr(:),':',naggr

  !
  ! COO_PROL should arrive here with local numbering
  !
  if (debug) write(0,*)  me,' ',trim(name),' Size check on entry New: ',&
       & coo_prol%get_fmt(),coo_prol%get_nrows(),coo_prol%get_ncols(),coo_prol%get_nzeros(),&
       & nrow,ntaggr,naggr

  call coo_prol%cp_to_fmt(csr_prol,info)

  if (debug) write(0,*) me,trim(name),' Product AxPROL ',&
       & a_csr%get_nrows(),a_csr%get_ncols(), csr_prol%get_nrows(), &
       & desc_a%get_local_rows(),desc_a%get_local_cols(),&
       & desc_ac%get_local_rows(),desc_a%get_local_cols()
  if (debug) flush(0)

  call psb_par_spspmm(a_csr,desc_a,csr_prol,acsr3,desc_ac,info)

  if (debug) write(0,*) me,trim(name),' Done AxPROL ',&
       & acsr3%get_nrows(),acsr3%get_ncols(), acsr3%get_nzeros(),&
       & desc_ac%get_local_rows(),desc_ac%get_local_cols()


  !
  ! Remember that RESTR must be built from PROL after halo extension,
  ! which is done above in psb_par_spspmm
  if (debug) write(0,*)me,' ',name,' No inp_restr, transposing prol ',&
       & csr_prol%get_nrows(),csr_prol%get_ncols(),csr_prol%get_nzeros()

  call csr_restr%cp_from_coo(coo_restr,info)

!!$      write(0,*)me,' ',name,' after transposition ',coo_restr%get_nrows(),coo_restr%get_ncols(),coo_restr%get_nzeros()

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv coo_restr')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  if (debug) write(0,*) me,trim(name),' Product RESTRxAP ',&
       & csr_restr%get_nrows(),csr_restr%get_ncols(), &
       & desc_ac%get_local_rows(),desc_a%get_local_cols(),&
       & acsr3%get_nrows(),acsr3%get_ncols()
  call psb_par_spspmm(csr_restr,desc_a,acsr3,ac_csr,desc_ac,info)
  call acsr3%free()


  call ac_csr%set_nrows(desc_ac%get_local_rows())
  call ac_csr%set_ncols(desc_ac%get_local_cols())
  call ac%mv_from(ac_csr)
  call ac%set_asb()

  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr
  ! write(0,*)  me,' ',trim(name),' Final AC newstyle ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros()

  call coo_prol%set_ncols(desc_ac%get_local_cols())
  if (debug) call check_coo(me,trim(name)//' Check 3 on coo_restr:',coo_restr)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done rap '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_c_rap

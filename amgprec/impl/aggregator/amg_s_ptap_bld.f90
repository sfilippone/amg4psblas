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
! File: amg_ptap_bld.F90
!
!  This routine does two things:
!   1. Computes the transpose of the prolongator Prol and stores
!      it explicitly into Restr
!   2. Computes  AC = P^T  A P    
!
!   Yes, doing two things at once is against received wisdom, but we are
!   trying to save time in the build procedure; maybe we will revisit this
!   decision. 
!
subroutine amg_s_ptap_bld(a_csr,desc_a,nlaggr,parms,ac,&
     & coo_prol,desc_ac,coo_restr,info,desc_ax)
  use psb_base_mod
  use amg_s_inner_mod
  use amg_s_base_aggregator_mod, amg_protect_name => amg_s_ptap_bld
  implicit none

  ! Arguments
  type(psb_s_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(inout)          :: desc_a
  integer(psb_lpk_), intent(inout)           :: nlaggr(:)
  type(amg_sml_parms), intent(inout)         :: parms 
  type(psb_s_coo_sparse_mat), intent(inout)  :: coo_prol, coo_restr
  type(psb_desc_type), intent(inout)         :: desc_ac
  type(psb_sspmat_type), intent(out)        :: ac
  integer(psb_ipk_), intent(out)             :: info
  type(psb_desc_type), intent(inout), optional :: desc_ax

  ! Local variables
  integer(psb_ipk_)  :: err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)  :: np, me, ndx
  character(len=40)  :: name
  type(psb_ls_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_s_csr_sparse_mat) :: acsr3, csr_prol, ac_csr, csr_restr
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nglob, ntaggr, naggrm1, naggrp1
  integer(psb_ipk_) :: nrow, ncol, nrl, nzl, ip, nzt, i, k
  integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza
  logical, parameter :: do_timings=.false., oldstyle=.false., debug=.false.  
  integer(psb_ipk_), save :: idx_spspmm=-1, idx_cpytrans1=-1, idx_cpytrans2=-1

  name='amg_ptap_bld'
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

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("PTAP_BLD: par_spspmm")
  if ((do_timings).and.(idx_cpytrans1==-1)) &
       & idx_cpytrans1 = psb_get_timer_idx("PTAP_BLD: cpy&trans1")
  if ((do_timings).and.(idx_cpytrans2==-1)) &
       & idx_cpytrans2 = psb_get_timer_idx("PTAP_BLD: cpy&trans2")

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

  if (do_timings) call psb_tic(idx_spspmm)
  call psb_par_spspmm(a_csr,desc_a,csr_prol,acsr3,desc_ac,info)
  if (do_timings) call psb_toc(idx_spspmm)  

  if (debug) write(0,*) me,trim(name),' Done AxPROL ',&
       & acsr3%get_nrows(),acsr3%get_ncols(), acsr3%get_nzeros(),&
       & desc_ac%get_local_rows(),desc_ac%get_local_cols()

  !
  ! Ok first product done.

  if (present(desc_ax)) then
    if (do_timings) call psb_tic(idx_cpytrans1)
    block 
      call coo_prol%cp_to_coo(coo_restr,info)
      call coo_restr%set_ncols(desc_ac%get_local_cols())
      call coo_restr%set_nrows(desc_a%get_local_rows())
      call psb_s_coo_glob_transpose(coo_restr,desc_a,info,desc_c=desc_ac,desc_rx=desc_ax)
      call coo_restr%set_nrows(desc_ac%get_local_rows())
      call coo_restr%set_ncols(desc_ax%get_local_cols())
    end block
    call csr_restr%cp_from_coo(coo_restr,info)
    if (do_timings) call psb_toc(idx_cpytrans1)
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
    if (do_timings) call psb_tic(idx_spspmm)      
    call psb_par_spspmm(csr_restr,desc_ax,acsr3,ac_csr,desc_ac,info)
    if (do_timings) call psb_toc(idx_spspmm)      
    call acsr3%free()

  else

    !
    ! Remember that RESTR must be built from PROL after halo extension,
    ! which is done above in psb_par_spspmm
    if (debug) write(0,*)me,' ',name,' No inp_restr, transposing prol ',&
         & csr_prol%get_nrows(),csr_prol%get_ncols(),csr_prol%get_nzeros()
    call csr_prol%mv_to_coo(coo_restr,info)
    if (debug) call check_coo(me,trim(name)//' Check 1 (before transp) on coo_restr:',coo_restr)

    call coo_restr%transp()
    nzl = coo_restr%get_nzeros()
    nrl = desc_ac%get_local_rows()
    call coo_restr%fix(info)
    i=coo_restr%get_nzeros()
    !
    ! Only keep local rows
    !
    search: do k=i,1,-1
      if (coo_restr%ia(k) <= nrl) then
        call coo_restr%set_nzeros(k)
        exit search
      end if
    end do search

    nzl  = coo_restr%get_nzeros()
    call coo_restr%set_nrows(desc_ac%get_local_rows())
    call coo_restr%set_ncols(desc_a%get_local_cols())
    if (debug) call check_coo(me,trim(name)//' Check 2 on coo_restr:',coo_restr)
    if (do_timings) call psb_tic(idx_cpytrans2)
        
    call csr_restr%cp_from_coo(coo_restr,info)

    if (do_timings) call psb_toc(idx_cpytrans2)
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
    if (do_timings) call psb_tic(idx_spspmm)      
    call psb_par_spspmm(csr_restr,desc_a,acsr3,ac_csr,desc_ac,info)
    if (do_timings) call psb_toc(idx_spspmm)      
    call acsr3%free()
  end if

  if (.not.desc_ac%is_asb()) call psb_cdasb(desc_ac,info)

  call ac_csr%set_nrows(desc_ac%get_local_rows())
  call ac_csr%set_ncols(desc_ac%get_local_cols())
  call ac%mv_from(ac_csr)
  call ac%set_asb()

  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr

  call coo_prol%set_ncols(desc_ac%get_local_cols())
  if (debug) call check_coo(me,trim(name)//' Check 3 on coo_restr:',coo_restr)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done ptap '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
contains
  subroutine check_coo(me,string,coo)
    implicit none
    integer(psb_ipk_) :: me
    type(psb_s_coo_sparse_mat) :: coo
    character(len=*) :: string
    integer(psb_lpk_) :: nr,nc,nz
    nr = coo%get_nrows()
    nc = coo%get_ncols()
    nz = coo%get_nzeros()
    write(0,*) me,string,nr,nc,&
         & minval(coo%ia(1:nz)),maxval(coo%ia(1:nz)),&
         & minval(coo%ja(1:nz)),maxval(coo%ja(1:nz))

  end subroutine check_coo
  
end subroutine amg_s_ptap_bld

subroutine amg_s_ls_ptap_bld(a_csr,desc_a,nlaggr,parms,ac,&
     & coo_prol,desc_ac,coo_restr,info,desc_ax)
  use psb_base_mod
  use amg_s_inner_mod
  use amg_s_base_aggregator_mod !, amg_protect_name => amg_s_ls_ptap_bld
  implicit none

  ! Arguments
  type(psb_s_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(inout)          :: desc_a
  integer(psb_lpk_), intent(inout)           :: nlaggr(:)
  type(amg_sml_parms), intent(inout)         :: parms 
  type(psb_ls_coo_sparse_mat), intent(inout)  :: coo_prol, coo_restr
  type(psb_desc_type), intent(inout)         :: desc_ac
  type(psb_lsspmat_type), intent(out)        :: ac
  integer(psb_ipk_), intent(out)             :: info
  type(psb_desc_type), intent(inout), optional :: desc_ax

  ! Local variables
  integer(psb_ipk_)   :: err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, ndx
  character(len=40)   :: name
  type(psb_ls_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_s_csr_sparse_mat) :: acsr3, csr_prol, ac_csr, csr_restr
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nrl, nzl, ip, &
       &  nzt, naggrm1, naggrp1, i, k
  integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza
  logical, parameter :: do_timings=.false., oldstyle=.false., debug=.false.  
  integer(psb_ipk_), save :: idx_spspmm=-1

  name='amg_ptap_bld'
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

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("SPMM_BLD: par_spspmm")

  naggr   = nlaggr(me+1)
  ntaggr  = sum(nlaggr)
  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1)) 

  !
  ! COO_PROL should arrive here with local numbering
  !
  if (debug) write(0,*)  me,' ',trim(name),' Size check on entry New: ',&
       & coo_prol%get_fmt(),coo_prol%get_nrows(),coo_prol%get_ncols(),coo_prol%get_nzeros(),&
       & nrow,ntaggr,naggr

  call coo_prol%cp_to_ifmt(csr_prol,info)

  if (debug) write(0,*) me,trim(name),' Product AxPROL ',&
       & a_csr%get_nrows(),a_csr%get_ncols(), csr_prol%get_nrows(), &
       & desc_a%get_local_rows(),desc_a%get_local_cols(),&
       & desc_ac%get_local_rows(),desc_a%get_local_cols()
  if (debug) flush(0)

  if (do_timings) call psb_tic(idx_spspmm)
  call psb_par_spspmm(a_csr,desc_a,csr_prol,acsr3,desc_ac,info)
  if (do_timings) call psb_toc(idx_spspmm)  

  if (debug) write(0,*) me,trim(name),' Done AxPROL ',&
       & acsr3%get_nrows(),acsr3%get_ncols(), acsr3%get_nzeros(),&
       & desc_ac%get_local_rows(),desc_ac%get_local_cols()

  !
  ! Ok first product done.

  if (present(desc_ax)) then
    block 
      type(psb_s_coo_sparse_mat) :: icoo_restr

      call coo_prol%cp_to_icoo(icoo_restr,info)
      call icoo_restr%set_ncols(desc_ac%get_local_cols())
      call icoo_restr%set_nrows(desc_a%get_local_rows())
      call psb_s_coo_glob_transpose(icoo_restr,desc_a,info,desc_c=desc_ac,desc_rx=desc_ax)
      call icoo_restr%set_nrows(desc_ac%get_local_rows())
      call icoo_restr%set_ncols(desc_ax%get_local_cols())
      call coo_restr%cp_from_icoo(icoo_restr,info)
    end block
    call csr_restr%cp_from_lcoo(coo_restr,info)


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
    if (do_timings) call psb_tic(idx_spspmm)      
    call psb_par_spspmm(csr_restr,desc_ax,acsr3,ac_csr,desc_ac,info)
    if (do_timings) call psb_toc(idx_spspmm)      
    call acsr3%free()

  else

    !
    ! Remember that RESTR must be built from PROL after halo extension,
    ! which is done above in psb_par_spspmm
    if (debug) write(0,*)me,' ',name,' No inp_restr, transposing prol ',&
         & csr_prol%get_nrows(),csr_prol%get_ncols(),csr_prol%get_nzeros()
    call csr_prol%mv_to_lcoo(coo_restr,info)
    if (debug) call check_coo(me,trim(name)//' Check 1 (before transp) on coo_restr:',coo_restr)

    call coo_restr%transp()
    nzl = coo_restr%get_nzeros()
    nrl = desc_ac%get_local_rows() 
    i=0
    !
    ! Only keep local rows
    !
    do k=1, nzl
      if ((1 <= coo_restr%ia(k)) .and.(coo_restr%ia(k) <= nrl)) then
        i = i+1
        coo_restr%val(i) = coo_restr%val(k)
        coo_restr%ia(i)  = coo_restr%ia(k)
        coo_restr%ja(i)  = coo_restr%ja(k)
      end if
    end do
    call coo_restr%set_nzeros(i)
    call coo_restr%fix(info) 
    nzl  = coo_restr%get_nzeros()
    call coo_restr%set_nrows(desc_ac%get_local_rows())
    call coo_restr%set_ncols(desc_a%get_local_cols())
    if (debug) call check_coo(me,trim(name)//' Check 2 on coo_restr:',coo_restr)
    call csr_restr%cp_from_lcoo(coo_restr,info)

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
    if (do_timings) call psb_tic(idx_spspmm)      
    call psb_par_spspmm(csr_restr,desc_a,acsr3,ac_csr,desc_ac,info)
    if (do_timings) call psb_toc(idx_spspmm)      
    call acsr3%free()
  end if

  if (.not.desc_ac%is_asb()) call psb_cdasb(desc_ac,info)

  call ac_csr%set_nrows(desc_ac%get_local_rows())
  call ac_csr%set_ncols(desc_ac%get_local_cols())
  call ac%mv_from(ac_csr)
  call ac%set_asb()

  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr

  call coo_prol%set_ncols(desc_ac%get_local_cols())
  if (debug) call check_coo(me,trim(name)//' Check 3 on coo_restr:',coo_restr)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done ptap '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
contains
  subroutine check_coo(me,string,coo)
    implicit none
    integer(psb_ipk_) :: me
    type(psb_ls_coo_sparse_mat) :: coo
    character(len=*) :: string
    integer(psb_lpk_) :: nr,nc,nz
    nr = coo%get_nrows()
    nc = coo%get_ncols()
    nz = coo%get_nzeros()
    write(0,*) me,string,nr,nc,&
         & minval(coo%ia(1:nz)),maxval(coo%ia(1:nz)),&
         & minval(coo%ja(1:nz)),maxval(coo%ja(1:nz))

  end subroutine check_coo
  
end subroutine amg_s_ls_ptap_bld

subroutine amg_ls_ptap_bld(a_csr,desc_a,nlaggr,parms,ac,&
     & coo_prol,desc_ac,coo_restr,info,desc_ax)
  use psb_base_mod
  use amg_s_inner_mod
  use amg_s_base_aggregator_mod!, amg_protect_name => amg_ls_ptap_bld
  implicit none

  ! Arguments
  type(psb_ls_csr_sparse_mat), intent(inout) :: a_csr
  type(psb_desc_type), intent(inout)           :: desc_a
  integer(psb_lpk_), intent(inout)           :: nlaggr(:)
  type(amg_sml_parms), intent(inout)         :: parms 
  type(psb_ls_coo_sparse_mat), intent(inout)  :: coo_prol, coo_restr
  type(psb_desc_type), intent(inout)         :: desc_ac
  type(psb_lsspmat_type), intent(out)        :: ac
  integer(psb_ipk_), intent(out)             :: info
  type(psb_desc_type), intent(inout), optional :: desc_ax

  ! Local variables
  integer(psb_ipk_)  :: err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)  :: np, me, ndx
  character(len=40)  :: name
  type(psb_ls_coo_sparse_mat) :: ac_coo, tmpcoo
  type(psb_ls_csr_sparse_mat) :: acsr3, csr_prol, ac_csr, csr_restr
  integer(psb_ipk_) :: debug_level, debug_unit, naggr
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nrl, nzl, ip, &
       &  nzt, naggrm1, naggrp1, i, k
  integer(psb_lpk_) ::  nrsave, ncsave, nzsave, nza
  logical, parameter :: do_timings=.true., oldstyle=.false., debug=.false.  
  integer(psb_ipk_), save :: idx_spspmm=-1

  name='amg_ptap_bld'
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

  if ((do_timings).and.(idx_spspmm==-1)) &
       & idx_spspmm = psb_get_timer_idx("SPMM_BLD: par_spspmm")

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

  if (do_timings) call psb_tic(idx_spspmm)
  call psb_par_spspmm(a_csr,desc_a,csr_prol,acsr3,desc_ac,info)
  if (do_timings) call psb_toc(idx_spspmm)  

  if (debug) write(0,*) me,trim(name),' Done AxPROL ',&
       & acsr3%get_nrows(),acsr3%get_ncols(), acsr3%get_nzeros(),&
       & desc_ac%get_local_rows(),desc_ac%get_local_cols()

  !
  ! Ok first product done.

  if (present(desc_ax)) then
    block 
      type(psb_s_coo_sparse_mat) :: icoo_restr

      call coo_prol%cp_to_icoo(icoo_restr,info)
      call icoo_restr%set_ncols(desc_ac%get_local_cols())
      call icoo_restr%set_nrows(desc_a%get_local_rows())
      call psb_s_coo_glob_transpose(icoo_restr,desc_a,info,desc_c=desc_ac,desc_rx=desc_ax)
      call icoo_restr%set_nrows(desc_ac%get_local_rows())
      call icoo_restr%set_ncols(desc_ax%get_local_cols())
      write(0,*) me,' ',trim(name),' check on glob_transpose 1: ',&
           & desc_a%get_local_cols(),desc_ax%get_local_cols(),icoo_restr%get_nzeros()
      if (desc_a%get_local_cols()<desc_ax%get_local_cols()) then
        write(0,*) me,' ',trim(name),' WARNING: GLOB_TRANSPOSE NEW INDICES '
      end if
      call coo_restr%cp_from_icoo(icoo_restr,info)
    end block
    call csr_restr%cp_from_coo(coo_restr,info)

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
    if (do_timings) call psb_tic(idx_spspmm)      
    call psb_par_spspmm(csr_restr,desc_ax,acsr3,ac_csr,desc_ac,info)
    if (do_timings) call psb_toc(idx_spspmm)      
    call acsr3%free()

  else

    !
    ! Remember that RESTR must be built from PROL after halo extension,
    ! which is done above in psb_par_spspmm
    if (debug) write(0,*)me,' ',name,' No inp_restr, transposing prol ',&
         & csr_prol%get_nrows(),csr_prol%get_ncols(),csr_prol%get_nzeros()
    call csr_prol%mv_to_coo(coo_restr,info)
    if (debug) call check_coo(me,trim(name)//' Check 1 (before transp) on coo_restr:',coo_restr)

    call coo_restr%transp()
    nzl = coo_restr%get_nzeros()
    nrl = desc_ac%get_local_rows() 
    i=0
    !
    ! Only keep local rows
    !
    do k=1, nzl
      if ((1 <= coo_restr%ia(k)) .and.(coo_restr%ia(k) <= nrl)) then
        i = i+1
        coo_restr%val(i) = coo_restr%val(k)
        coo_restr%ia(i)  = coo_restr%ia(k)
        coo_restr%ja(i)  = coo_restr%ja(k)
      end if
    end do
    call coo_restr%set_nzeros(i)
    call coo_restr%fix(info) 
    nzl  = coo_restr%get_nzeros()
    call coo_restr%set_nrows(desc_ac%get_local_rows())
    call coo_restr%set_ncols(desc_a%get_local_cols())
    if (debug) call check_coo(me,trim(name)//' Check 2 on coo_restr:',coo_restr)
    call csr_restr%cp_from_coo(coo_restr,info)

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
    if (do_timings) call psb_tic(idx_spspmm)      
    call psb_par_spspmm(csr_restr,desc_a,acsr3,ac_csr,desc_ac,info)
    if (do_timings) call psb_toc(idx_spspmm)      
    call acsr3%free()
  end if
  
  if (.not.desc_ac%is_asb()) call psb_cdasb(desc_ac,info)

  call ac_csr%set_nrows(desc_ac%get_local_rows())
  call ac_csr%set_ncols(desc_ac%get_local_cols())
  call ac%mv_from(ac_csr)
  call ac%set_asb()

  if (debug) write(0,*)  me,' ',trim(name),' After  mv_from',psb_get_errstatus()
  if (debug) write(0,*)  me,' ',trim(name),' ',ac%get_fmt(),ac%get_nrows(),ac%get_ncols(),ac%get_nzeros(),naggr,ntaggr

  call coo_prol%set_ncols(desc_ac%get_local_cols())
  if (debug) call check_coo(me,trim(name)//' Check 3 on coo_restr:',coo_restr)


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done ptap '

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains
  subroutine check_coo(me,string,coo)
    implicit none
    integer(psb_ipk_) :: me
    type(psb_ls_coo_sparse_mat) :: coo
    character(len=*) :: string
    integer(psb_lpk_) :: nr,nc,nz
    nr = coo%get_nrows()
    nc = coo%get_ncols()
    nz = coo%get_nzeros()
    write(0,*) me,string,nr,nc,&
         & minval(coo%ia(1:nz)),maxval(coo%ia(1:nz)),&
         & minval(coo%ja(1:nz)),maxval(coo%ja(1:nz))

  end subroutine check_coo
  
end subroutine amg_ls_ptap_bld

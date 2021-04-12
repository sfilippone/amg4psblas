!   !
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
!  moved here from
!
!                             AMG4PSBLAS  Extensions
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
! File: amg_s_parmatch_aggregator_tprol.f90
!
! Subroutine: amg_s_parmatch_aggregator_tprol
! Version:    real
!
!

subroutine  amg_s_parmatch_aggregator_build_tprol(ag,parms,ag_data,&
     & a,desc_a,ilaggr,nlaggr,t_prol,info)
  use psb_base_mod
  use amg_base_prec_type
  use amg_s_inner_mod
  use amg_s_parmatch_aggregator_mod, amg_protect_name => amg_s_parmatch_aggregator_build_tprol
  use iso_c_binding
  implicit none
  class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
  type(amg_sml_parms), intent(inout)   :: parms
  type(amg_saggr_data), intent(in)     :: ag_data
  type(psb_sspmat_type), intent(inout) :: a
  type(psb_desc_type), intent(inout)   :: desc_a
  integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
  type(psb_lsspmat_type), intent(out)  :: t_prol
  integer(psb_ipk_), intent(out)      :: info


  ! Local variables
  real(psb_spk_), allocatable    :: tmpw(:), tmpwnxt(:)
  integer(psb_lpk_), allocatable :: ixaggr(:), nxaggr(:), tlaggr(:), ivr(:)
  type(psb_sspmat_type)          :: a_tmp
  integer(c_int) :: match_algorithm, n_sweeps, max_csize, max_nlevels
  character(len=40)    :: name, ch_err
  character(len=80)    :: fname, prefix_
  type(psb_ctxt_type)  :: ictxt
  integer(psb_ipk_)    :: np, me
  integer(psb_ipk_)    :: err_act, ierr
  integer(psb_ipk_)    :: debug_level, debug_unit
  integer(psb_ipk_)    :: i, j, k, nr, nc
  integer(psb_lpk_)    :: isz, num_pcols, nrac, ncac, lname, nz, x_sweeps, csz
  integer(psb_lpk_)    :: psz, sizes(4)
  type(psb_s_csr_sparse_mat), target  :: csr_prol, csr_pvi, csr_prod_res, acsr
  type(psb_ls_csr_sparse_mat), target :: lcsr_prol
  type(psb_desc_type), allocatable    :: desc_acv(:)
  type(psb_ls_coo_sparse_mat)         :: tmpcoo, transp_coo
  type(psb_sspmat_type), allocatable  :: acv(:)
  type(psb_sspmat_type), allocatable  :: prolv(:), restrv(:)
  type(psb_lsspmat_type)              :: tmp_prol, tmp_pg, tmp_restr
  type(psb_desc_type)                 :: tmp_desc_ac, tmp_desc_ax, tmp_desc_p
  integer(psb_ipk_), save             :: idx_mboxp=-1, idx_spmmbld=-1, idx_sweeps_mult=-1
  logical, parameter :: dump=.false., do_timings=.true., debug=.false., &
       & dump_prol_restr=.false.

  name='s_parmatch_tprol'
  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)
  if (psb_get_errstatus().ne.0) then
    write(0,*) me,trim(name),' Err_status :',psb_get_errstatus()
    return
  end if
  if (debug) write(0,*) me,trim(name),' Start '
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info  = psb_success_

  if ((do_timings).and.(idx_mboxp==-1))       &
       & idx_mboxp = psb_get_timer_idx("PMC_TPROL: MatchBoxP")
  if ((do_timings).and.(idx_spmmbld==-1))     &
       & idx_spmmbld = psb_get_timer_idx("PMC_TPROL: spmm_bld")
  if ((do_timings).and.(idx_sweeps_mult==-1)) &
       & idx_sweeps_mult = psb_get_timer_idx("PMC_TPROL: sweeps_mult")


  call amg_check_def(parms%ml_cycle,'Multilevel cycle',&
       &   amg_mult_ml_,is_legal_ml_cycle)
  call amg_check_def(parms%par_aggr_alg,'Aggregation',&
       &   amg_coupled_aggr_,is_legal_coupled_par_aggr_alg)
  call amg_check_def(parms%aggr_ord,'Ordering',&
       &   amg_aggr_ord_nat_,is_legal_ml_aggr_ord)
  call amg_check_def(parms%aggr_thresh,'Aggr_Thresh',szero,is_legal_s_aggr_thrs)
  match_algorithm = ag%matching_alg
  n_sweeps        = ag%n_sweeps
  if (2**n_sweeps /= ag%orig_aggr_size) then
    if (me == 0) then
      write(debug_unit, *) 'Warning: AGGR_SIZE reset to value ',2**n_sweeps
    end if
  end if
  if (ag%max_csize > 0) then
    max_csize       = ag%max_csize
  else
    max_csize       = ag_data%min_coarse_size
  end if
  if (ag%max_nlevels > 0) then
    max_nlevels     = ag%max_nlevels
  else
    max_nlevels = ag_data%max_levs
  end if
  if (.true.) then
    block
      integer(psb_ipk_) :: ipv(2)
      ipv(1) = max_csize
      ipv(2) = n_sweeps
      call psb_bcast(ictxt,ipv)
      max_csize = ipv(1)
      n_sweeps  = ipv(2)
    end block
  else
    call psb_bcast(ictxt,max_csize)
    call psb_bcast(ictxt,n_sweeps)
  end if
  if (n_sweeps /= ag%n_sweeps) then
    write(0,*) me,' Inconsistent N_SWEEPS ',n_sweeps,ag%n_sweeps
  end if
!!$  if (me==0) write(0,*) 'Matching sweeps: ',n_sweeps
  n_sweeps = max(1,n_sweeps)
  if (debug) write(0,*) me,' Copies, with n_sweeps: ',n_sweeps,max_csize
  if (ag%unsmoothed_hierarchy.and.allocated(ag%base_a)) then
    call ag%base_a%cp_to(acsr)
    if (ag%do_clean_zeros) call acsr%clean_zeros(info)
    nr = acsr%get_nrows()
    if (psb_size(ag%w) < nr) call ag%bld_default_w(nr)
    isz = acsr%get_ncols()

    call psb_realloc(isz,ixaggr,info)
    if (info == psb_success_) &
         & allocate(acv(0:n_sweeps), desc_acv(0:n_sweeps),&
         &  prolv(n_sweeps), restrv(n_sweeps),stat=info)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999

    end if


    call acv(0)%mv_from(acsr)
    call ag%base_desc%clone(desc_acv(0),info)

  else
    call a%cp_to(acsr)
    if (ag%do_clean_zeros) call acsr%clean_zeros(info)
    nr = acsr%get_nrows()
    if (psb_size(ag%w) < nr) call ag%bld_default_w(nr)
    isz = acsr%get_ncols()

    call psb_realloc(isz,ixaggr,info)
    if (info == psb_success_) &
         & allocate(acv(0:n_sweeps), desc_acv(0:n_sweeps),&
         &  prolv(n_sweeps), restrv(n_sweeps),stat=info)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999

    end if


    call acv(0)%mv_from(acsr)
    call desc_a%clone(desc_acv(0),info)
  end if

  nrac = desc_acv(0)%get_local_rows()
  ncac = desc_acv(0)%get_local_cols()
  if (debug) write(0,*) me,' On input to level: ',nrac, ncac
  if (allocated(ag%prol)) then
    call ag%prol%free()
    deallocate(ag%prol)
  end if
  if (allocated(ag%restr)) then
    call ag%restr%free()
    deallocate(ag%restr)
  end if

  if (dump) then
    block
      type(psb_lsspmat_type) :: lac
      ivr = desc_acv(0)%get_global_indices(owned=.false.)
      prefix_ = "input_a"
      lname = len_trim(prefix_)
      fname = trim(prefix_)
      write(fname(lname+1:lname+9),'(a,i3.3,a)') '_p',me, '.mtx'
      call acv(0)%print(fname,head='Debug aggregates')
      call lac%cp_from(acv(0))
      write(fname(lname+1:lname+13),'(a,i3.3,a)') '_p',me, '-glb.mtx'
      call lac%print(fname,head='Debug aggregates',iv=ivr)
      call lac%free()
    end block
  end if

  call psb_geall(tmpw,desc_acv(0),info)

  tmpw(1:nr) = ag%w(1:nr)

  call psb_geasb(tmpw,desc_acv(0),info)

  if (debug) then
    call psb_barrier(ictxt)
    if (me == 0) write(0,*) 'N_sweeps ',n_sweeps,nr,desc_acv(0)%is_ok(),max_csize
  end if

  !
  ! Prepare ag%ac, ag%desc_ac, ag%prol, ag%restr to enable
  !         shortcuts in mat_bld and mat_asb
  ! and ag%desc_ax which will be needed in backfix.
  !
  x_sweeps = -1
  sweeps_loop: do i=1, n_sweeps
    if (debug) then
      call psb_barrier(ictxt)
      if (me==0) write(0,*) me,trim(name),' Start sweeps_loop iteration:',i,' of ',n_sweeps
    end if

    !
    ! Building prol and restr because this algorithm is not decoupled
    ! On exit from matchbox_build_prol,  prolv(i) is in global numbering
    !
    !
    if (debug) write(0,*) me,' Into matchbox_build_prol ',info
    if (do_timings) call psb_tic(idx_mboxp)
    call smatchboxp_build_prol(tmpw,acv(i-1),desc_acv(i-1),ixaggr,nxaggr,tmp_prol,info,&
         & symmetrize=ag%need_symmetrize,reproducible=ag%reproducible_matching)
    if (do_timings) call psb_toc(idx_mboxp)
    if (debug) write(0,*) me,' Out from matchbox_build_prol ',info
    if (psb_errstatus_fatal())  write(0,*)me,trim(name),'Error fatal on exit bld_tprol',info


    if (debug) then
      call psb_barrier(ictxt)
!!$      write(0,*) name,' Call spmm_bld sweep:',i,n_sweeps
      if (me==0) write(0,*) me,trim(name),' Calling spmm_bld  NSW>1:',i,&
           & desc_acv(i-1)%get_local_rows(),desc_acv(i-1)%get_local_cols(),&
           & desc_acv(i-1)%get_global_rows()
    end if
    if (i == n_sweeps) call tmp_prol%clone(tmp_pg,info)
    if (do_timings) call psb_tic(idx_spmmbld)
    !
    !  On entry, prolv(i) is in global numbering,
    !
    call amg_s_parmatch_spmm_bld_ov(acv(i-1),desc_acv(i-1),ixaggr,nxaggr,parms,&
         & acv(i),desc_acv(i), prolv(i),restrv(1),tmp_prol,info)
    if (psb_errstatus_fatal())  write(0,*)me,trim(name),'Error fatal on exit from bld_ov(i)',info
    if (debug) then
      call psb_barrier(ictxt)
      if (me==0) write(0,*) me,trim(name),' Done spmm_bld:',i
    end if

    if (do_timings) call psb_toc(idx_spmmbld)
    ! Keep a copy of prolv(i) in global numbering for the time being, will
    ! need it to build the final
    ! if (i == n_sweeps) call prolv(i)%clone(tmp_prol,info)
    call ag%inner_mat_asb(parms,acv(i-1),desc_acv(i-1),&
         & acv(i),desc_acv(i),prolv(i),restrv(1),info)

    if (debug) then
      call psb_barrier(ictxt)
      if (me==0) write(0,*) me,trim(name),' Done mat_asb:',i,sum(nxaggr),max_csize,info
      csz = sum(nxaggr)
      call psb_bcast(ictxt,csz)
      if (csz /= sum(nxaggr)) write(0,*) me,trim(name),' Mismatch matasb',&
           & csz,sum(nxaggr),max_csize
    end if
    if (psb_errstatus_fatal())  write(0,*)me,trim(name),'Error fatal on entry to tmpwnxt 2'


    !
    ! Fix wnxt
    !
    if (info == 0) call psb_geall(tmpwnxt,desc_acv(i),info)
    if (info == 0) call psb_geasb(tmpwnxt,desc_acv(i),info,scratch=.true.)
    if (info == 0) call psb_halo(tmpw,desc_acv(i-1),info)
!!$      write(0,*) trestr%get_nrows(),size(tmpwnxt),trestr%get_ncols(),size(tmpw)

    if (info == 0) call psb_csmm(sone,restrv(1),tmpw,szero,tmpwnxt,info)

    if (info /= psb_success_) then
      write(0,*)me,trim(name),'Error from mat_asb/tmpw ',info
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='mat_asb 2')
      goto 9999
    end if


    if (i == 1) then
      nrac = desc_acv(1)%get_local_rows()
!!$      write(0,*) 'Copying output w_nxt ',nrac
      call psb_realloc(nrac,ag%w_nxt,info)
      ag%w_nxt(1:nrac) = tmpwnxt(1:nrac)
      !
      !  ILAGGR is fixed later on, but
      !  get a copy in case of an early exit
      !
      call psb_safe_ab_cpy(ixaggr,ilaggr,info)
    end if
    call psb_safe_ab_cpy(nxaggr,nlaggr,info)
    call move_alloc(tmpwnxt,tmpw)
    if (debug) then
      if (csz /= sum(nlaggr)) write(0,*) me,trim(name),' Mismatch 2 matasb',&
           & csz,sum(nlaggr),max_csize, info
    end if
    call acv(i-1)%free()
    if ((sum(nlaggr) <= max_csize).or.(any(nlaggr==0))) then
      x_sweeps = i
      exit sweeps_loop
    end if
    if (debug) then
      call psb_barrier(ictxt)
      if (me==0) write(0,*) me,trim(name),' Done sweeps_loop iteration:',i,' of ',n_sweeps
    end if

  end do sweeps_loop

  if (debug) then
    call psb_barrier(ictxt)
    if (me==0) write(0,*) me,trim(name),' Done sweeps_loop:',x_sweeps
  end if
  if (x_sweeps<=0) x_sweeps = n_sweeps

  if (do_timings) call psb_tic(idx_sweeps_mult)
  !
  ! Ok, now we have all the prolongators, including the last one in global numbering.
  ! Build the product of all prolongators. Need a tmp_desc_ax
  ! which is correct but most of the time overdimensioned
  !
  if (.not.allocated(ag%desc_ax)) allocate(ag%desc_ax)
  !
  block
    integer(psb_ipk_) :: i, nnz
    integer(psb_lpk_) :: ncol, ncsave
    if (.not.allocated(ag%ac)) allocate(ag%ac)
    if (.not.allocated(ag%desc_ac)) allocate(ag%desc_ac)
    call desc_acv(x_sweeps)%clone(ag%desc_ac,info)
    call desc_acv(x_sweeps)%free(info)
    call acv(x_sweeps)%move_alloc(ag%ac,info)
    if (.not.allocated(ag%prol)) allocate(ag%prol)
    if (.not.allocated(ag%restr)) allocate(ag%restr)

    call psb_cd_reinit(ag%desc_ac,info)
    ncsave =  ag%desc_ac%get_global_rows()
    !
    ! Note: prolv(i) is already in local numbering
    ! because of the call to mat_asb in the loop above.
    !
    call prolv(x_sweeps)%mv_to(csr_prol)
    if (debug) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*) 'Enter prolongator product loop ',x_sweeps
    end if

    do i=x_sweeps-1, 1, -1
      call prolv(i)%mv_to(csr_pvi)
      if (psb_errstatus_fatal()) write(0,*) me,' Fatal error in prolongator loop 1'
      call psb_par_spspmm(csr_pvi,desc_acv(i),csr_prol,csr_prod_res,ag%desc_ac,info)
      if ((info /=0).or.psb_errstatus_fatal()) write(0,*) me,' Fatal error in prolongator loop 2',info
      call csr_pvi%free()
      call csr_prod_res%mv_to_fmt(csr_prol,info)
      if ((info /=0).or.psb_errstatus_fatal()) write(0,*) me,' Fatal error in prolongator loop 3',info
      call csr_prol%set_ncols(ag%desc_ac%get_local_cols())
      if ((info /=0).or.psb_errstatus_fatal()) write(0,*) me,' Fatal error in prolongator loop 4'
    end do
    call csr_prol%mv_to_lfmt(lcsr_prol,info)
    nnz = lcsr_prol%get_nzeros()
    call ag%desc_ac%l2gip(lcsr_prol%ja(1:nnz),info)
    call lcsr_prol%set_ncols(ncsave)
    if (debug) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*) 'Done prolongator product loop ',x_sweeps
    end if
    !
    ! Fix ILAGGR here by copying from CSR_PROL%JA
    !
    block
      integer(psb_ipk_) :: nr
      nr = lcsr_prol%get_nrows()
      if (nnz /= nr) then
        write(0,*) me,name,' Issue with prolongator? ',nr,nnz
      end if
      call psb_realloc(nr,ilaggr,info)
      ilaggr(1:nnz) = lcsr_prol%ja(1:nnz)
    end block
    call tmp_prol%mv_from(lcsr_prol)
    call psb_cdasb(ag%desc_ac,info)
    call ag%ac%set_ncols(ag%desc_ac%get_local_cols())
  end block

  call tmp_prol%move_alloc(t_prol,info)
  call t_prol%set_ncols(ag%desc_ac%get_local_cols())
  call t_prol%set_nrows(desc_acv(0)%get_local_rows())

  nrac = ag%desc_ac%get_local_rows()
  ncac = ag%desc_ac%get_local_cols()
  call psb_realloc(nrac,ag%w_nxt,info)
  ag%w_nxt(1:nrac) = tmpw(1:nrac)


  if (do_timings) call psb_toc(idx_sweeps_mult)

  if (debug) then
    call psb_barrier(ictxt)
    if (me == 0) write(0,*) 'Out of build loop ',x_sweeps,': Output size:',sum(nlaggr)
  end if


  !call psb_set_debug_level(0)
  if (dump) then
    block
      ivr = desc_acv(x_sweeps)%get_global_indices(owned=.false.)
      prefix_ = "final_ac"
      lname = len_trim(prefix_)
      fname = trim(prefix_)
      write(fname(lname+1:lname+9),'(a,i3.3,a)') '_p',me, '.mtx'
      call acv(x_sweeps)%print(fname,head='Debug aggregates')
      write(fname(lname+1:lname+13),'(a,i3.3,a)') '_p',me, '-glb.mtx'
      call acv(x_sweeps)%print(fname,head='Debug aggregates',iv=ivr)
      prefix_ = "final_tp"
      lname = len_trim(prefix_)
      fname = trim(prefix_)
      write(fname(lname+1:lname+9),'(a,i3.3,a)') '_p',me, '.mtx'
      call t_prol%print(fname,head='Tentative prolongator')
    end block
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_bootCMatch_if')
    goto 9999
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
contains
  subroutine do_l1_jacobi(nsweeps,w,a,desc_a)
    integer(psb_ipk_), intent(in) :: nsweeps
    real(psb_dpk_), intent(inout) :: w(:)
    type(psb_sspmat_type), intent(in) :: a
    type(psb_desc_type), intent(in)   :: desc_a

  end subroutine do_l1_jacobi
end subroutine amg_s_parmatch_aggregator_build_tprol

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
!         software without specific prior written permission.
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
! File: amg_d_hierarchy_bld.f90
!
! Subroutine: amg_d_hierarchy_bld
! Version:    real
!
!  This routine builds the preconditioner according to the requirements made by
!  the user trough the subroutines amg_precinit and amg_precset.
!  
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level,
!  (for more details see the description of amg_Tonelev_type in amg_prec_type.f90).
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level. No transfer operators are associated to level 1.
! 
!
! Arguments:
!    a       -  type(psb_dspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_dprec_type), input/output.
!               The preconditioner data structure; upon exit it contains 
!               the multilevel hierarchy of prolongators, restrictors
!               and coarse matrices.
!    info    -  integer, output.
!               Error code.              
!  
subroutine amg_d_hierarchy_bld(a,desc_a,prec,info,cpymat)

  use psb_base_mod
  use amg_d_inner_mod
  use amg_d_prec_mod, amg_protect_name => amg_d_hierarchy_bld

  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(inout), target         :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  class(amg_dprec_type),intent(inout),target           :: prec
  integer(psb_ipk_), intent(out)                       :: info
  logical, intent(in), optional                        :: cpymat

  ! Local Variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me,np
  integer(psb_ipk_)   :: err,i,k, err_act, iszv, newsz,&
       & nplevs, mxplevs, level
  integer(psb_lpk_) :: iaggsize, casize, mncsize, mncszpp
  real(psb_dpk_)     :: mnaggratio, sizeratio, athresh, aomega
  class(amg_d_base_smoother_type), allocatable :: coarse_sm, med_sm, &
       & med_sm2, coarse_sm2
  class(amg_d_base_aggregator_type), allocatable :: tmp_aggr
  type(amg_dml_parms)              :: medparms, coarseparms
  integer(psb_lpk_), allocatable   :: ilaggr(:), nlaggr(:)
  type(psb_ldspmat_type)            :: op_prol
  type(amg_d_onelev_type), allocatable :: tprecv(:)
  logical :: cpymat_
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name
  character(len=40)  :: ch_err
  integer(psb_ipk_), save  :: idx_bldtp=-1, idx_matasb=-1
  logical, parameter :: do_timings=.false.
  logical             :: stop_hierarchy_loop
  type(psb_ctxt_type) :: lctxt
  integer(psb_ipk_)   :: lme,lnp

  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'amg_d_hierarchy_bld'
  info = psb_success_
  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  prec%ctxt = ctxt
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  if ((do_timings).and.(idx_bldtp==-1))       &
       & idx_bldtp = psb_get_timer_idx("BLD_HIER: bld_tprol")
  if ((do_timings).and.(idx_matasb==-1))     &
       & idx_matasb = psb_get_timer_idx("BLD_HIER: mmat_asb")
  !

  if (.not.allocated(prec%precv)) then 
    !! Error: should have called amg_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if
  cpymat_ = .false.
  if (present(cpymat)) cpymat_ = cpymat

  !
  ! Check to ensure all procs have the same 
  !   
  newsz      = -1
  mxplevs    = prec%ag_data%max_levs
  mnaggratio = prec%ag_data%min_cr_ratio
  mncsize    = prec%ag_data%min_coarse_size
  mncszpp    = prec%ag_data%min_coarse_size_per_process
  iszv       = prec%get_nlevs()
  call psb_bcast(ctxt,iszv)
  call psb_bcast(ctxt,mncsize)
  call psb_bcast(ctxt,mncszpp)
  call psb_bcast(ctxt,mxplevs)
  call psb_bcast(ctxt,mnaggratio)
  if (mncsize /= prec%ag_data%min_coarse_size) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_coarse_size')
    goto 9999
  end if
  if (mncszpp /= prec%ag_data%min_coarse_size_per_process) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_coarse_size_per_process')
    goto 9999
  end if
  if (mxplevs /= prec%ag_data%max_levs) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent max_levs')
    goto 9999
  end if
  if (mnaggratio /= prec%ag_data%min_cr_ratio) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_cr_ratio')
    goto 9999
  end if
  if (iszv /= prec%get_nlevs()) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if

  if (iszv < 1) then
    !
    ! This is wrong, cannot be size <1
    !
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  if (iszv == 1) then
    !
    ! This is OK, since it may be called by the user even if there
    ! is only one level
    !
    if (cpymat_) then
      call a%clone(prec%precv(1)%ac,info)
      call desc_a%clone(prec%precv(1)%desc_ac,info)
      prec%precv(1)%base_a    => prec%precv(1)%ac
      prec%precv(1)%base_desc => prec%precv(1)%desc_ac
    else
      prec%precv(1)%base_a    => a
      prec%precv(1)%base_desc => desc_a
    end if

    call psb_erractionrestore(err_act)
    return
  endif

  !
  ! The strategy:
  ! 1. The maximum number of levels should be already encoded in the
  !    size of the array;
  ! 2. If the user did not specify anything, then a default coarse size
  !    is generated, and the number of levels is set to the maximum;
  ! 3. If the size of the array is different from target number of levels,
  !    reallocate;
  ! 4. Build the matrix hierarchy, stopping early if either the target
  !    coarse size is hit, or the gain falls below the min_cr_ratio
  !    threshold.
  !
  if ((mncszpp < 0).and.(mncsize<0)) then 
    mncszpp = 200
    prec%ag_data%min_coarse_size_per_process = mncszpp
  end if
  if (mncszpp > 0) then
    casize = mncszpp*np
    if (casize > huge(ione)) casize = huge(ione)
  else
    if (mncsize < np) then
      if (me == 0) write(0,*) &
           & 'Warning: resetting coarse size to NP (1 variable per process)'
      mncsize = np
    end if
    casize = mncsize
  end if
  prec%ag_data%target_coarse_size = casize
  nplevs = max(itwo,mxplevs)

  !
  ! The coarse parameters will be needed later
  !
  coarseparms = prec%precv(iszv)%parms
  call save_smoothers(prec%precv(iszv),coarse_sm,coarse_sm2,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,a_err='Base level precbuild.')
    goto 9999
  end if
  !
  ! First set desired number of levels if different from default. 
  !
  if (iszv /= nplevs) then
    allocate(tprecv(nplevs),stat=info)
    ! First all existing levels
    do i=1, min(iszv,nplevs) - 1
      if (info == 0) tprecv(i)%parms = prec%precv(i)%parms
      if (info == 0) call restore_smoothers(tprecv(i),&
           & prec%precv(i)%sm,prec%precv(i)%sm2a,info)
      if (info == 0) call move_alloc(prec%precv(i)%aggr,tprecv(i)%aggr)
    end do
    if (iszv < nplevs) then 
      ! Further intermediates, if needed
      allocate(tmp_aggr,source=tprecv(iszv-1)%aggr,stat=info)
      medparms = prec%precv(iszv-1)%parms
      call save_smoothers(prec%precv(iszv-1),med_sm,med_sm2,info)
      do i=iszv, nplevs - 1
        if (info == 0) tprecv(i)%parms = medparms
        if (info == 0) call restore_smoothers(tprecv(i),med_sm,med_sm2,info)
        if ((info == 0).and..not.allocated(tprecv(i)%aggr))&
             & allocate(tprecv(i)%aggr,source=tmp_aggr,stat=info)
      end do
      deallocate(tmp_aggr,stat=info)
    end if

    ! Then coarse
    if (info == 0) tprecv(nplevs)%parms = coarseparms
    if (info == 0) call restore_smoothers(tprecv(nplevs),coarse_sm,coarse_sm2,info)
    if (info == 0) then
      if (nplevs <= iszv) then 
        allocate(tprecv(nplevs)%aggr,source=prec%precv(nplevs)%aggr,stat=info)
      else
        allocate(tmp_aggr,source=tprecv(nplevs-1)%aggr,stat=info)
        call move_alloc(tmp_aggr,tprecv(nplevs)%aggr)
      end if
    end if
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif

    do i=1,iszv
      call prec%precv(i)%free(info)
    end do
    call move_alloc(tprecv,prec%precv)
    call prec%set_nlevs(nplevs)
    iszv = prec%get_nlevs()
  end if

  !
  ! Finest level first; create a GEN_BLOCK
  ! copy of the descriptor. 
  !
  if (cpymat_) then
    call a%clone(prec%precv(1)%ac,info)
    prec%precv(1)%base_a    => prec%precv(1)%ac
  else
    prec%precv(1)%base_a    => a
  end if
  call psb_cd_renum_block(desc_a,prec%precv(1)%desc_ac,info)
  prec%precv(1)%base_desc => prec%precv(1)%desc_ac


  !
  ! Main build loop
  !   
  newsz = 0
  stop_hierarchy_loop = .false. 
  array_build_loop: do i=2, iszv
    !
    ! Check on the iprcparm contents: they should be the same
    ! on all processes.
    !
    call psb_bcast(ctxt,prec%precv(i)%parms)
    !
    ! Get current context: might have performed remapping
    ! 
    lctxt = prec%precv(i-1)%base_desc%get_ctxt()
    call psb_info(lctxt,lme,lnp)
!!$    write(0,*) 'Check at level',i,lme,lnp
    !
    ! Sanity checks on the parameters
    !
    if (i<iszv) then 
      !
      ! A replicated matrix only makes sense at the coarsest level
      !
      call amg_check_def(prec%precv(i)%parms%coarse_mat,'Coarse matrix',&
           &   amg_distr_mat_,is_distr_ml_coarse_mat)
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Calling mlprcbld at level  ',i
    !
    ! Build the tentative mapping between levels i-1 and i
    ! and the matrixat level i
    !
    if (do_timings) call psb_tic(idx_bldtp)
    if (info == psb_success_)&
         & call prec%precv(i)%bld_tprol(prec%precv(i-1)%base_a,&
         & prec%precv(i-1)%base_desc,&
         & ilaggr,nlaggr,op_prol,prec%ag_data,info)
    if (do_timings) call psb_toc(idx_bldtp)
    if (info /= psb_success_) then
      write(ch_err,'(a,i7)') 'Map bld fail @ level ',i
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err=ch_err)
      goto 9999
    endif

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Return from ',i,' call to bld_tprol', info      
    !
    ! Save op_prol just in case
    !
    call op_prol%clone(prec%precv(i)%tprol,info)

    !
    ! Check for early termination of aggregation loop. 
    !
    if (i == 2) then 
      call amg_d_hierarchy_bld_cmp_newsz(i,iszv,&
           & desc_a%get_global_rows(),&
           & nlaggr,casize,mnaggratio,sizeratio,newsz)
    else
      call amg_d_hierarchy_bld_cmp_newsz(i,iszv,&
           & sum(prec%precv(i-1)%linmap%naggr),&
           & nlaggr,casize,mnaggratio,sizeratio,newsz)
    end if
    prec%precv(i)%szratio = sizeratio
    call psb_bcast(ctxt,newsz)

    !
    ! Handle reallocation, if needed, and then mat_asb to polish off the
    ! construction
    ! 
    if (newsz > 0) then
      !
      ! This is awkward, we are saving the aggregation parms, for the sake
      ! of distr/repl matrix at coarse level. Should be rethought.
      !
      athresh =  prec%precv(newsz)%parms%aggr_thresh
      aomega  =  prec%precv(newsz)%parms%aggr_omega_val
      if (info == 0) prec%precv(newsz)%parms = coarseparms
      prec%precv(newsz)%parms%aggr_thresh    =  athresh
      prec%precv(newsz)%parms%aggr_omega_val =  aomega 

      if (info == 0) call restore_smoothers(prec%precv(newsz),&
           & coarse_sm,coarse_sm2,info)
      if (newsz < i) then
        !
        ! We are going back and revisit a previous level;
        ! recover the aggregation.
        !
        ilaggr = prec%precv(newsz)%linmap%iaggr
        nlaggr = prec%precv(newsz)%linmap%naggr
        call prec%precv(newsz)%tprol%clone(op_prol,info)
      end if
      if (do_timings) call psb_tic(idx_matasb)      
      if (info == psb_success_) call prec%precv(newsz)%mat_asb( &
           & prec%precv(newsz-1)%base_a,prec%precv(newsz-1)%base_desc,&
           & ilaggr,nlaggr,op_prol,info)
      if (do_timings) call psb_toc(idx_matasb)      
      if (info /= 0) then 
        write(ch_err,'(a,i7)') 'Mat asb fail @ level ',newsz
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err=ch_err)
        goto 9999
      endif
!!$      write(0,*) ' Early exit  of array_build_loop',i,iszv,info,&
      level = newsz
      stop_hierarchy_loop = .true. 
    else
      if (do_timings) call psb_tic(idx_matasb)
      if (info == psb_success_) call prec%precv(i)%mat_asb(&
           & prec%precv(i-1)%base_a,prec%precv(i-1)%base_desc,&
           & ilaggr,nlaggr,op_prol,info)
      if (do_timings) call psb_toc(idx_matasb)      
      level = i
    end if

    !
    ! Do we want to remap onto a smaller subset of processes?
    ! Will need a more sophisticated policy
    !
    block
      type(psb_ctxt_type) :: lctxt
      integer(psb_ipk_)   :: lme,lnp
      lctxt = prec%precv(level)%desc_ac%get_ctxt()
      call psb_info(lctxt,lme,lnp)
      if (amg_d_policy_do_remap(lctxt,level,sum(nlaggr))) then 
!!$        write(0,*) ' Context on remapping ',lme,lnp
        if ((lme >=0).and.(lnp>=2)) then 
          associate(lv=>prec%precv(level), rmp => prec%precv(level)%remap_data)
            call lv%desc_ac%clone(rmp%desc_ac_pre_remap,info)
            call lv%ac%clone(rmp%ac_pre_remap,info)
!!$            write(0,*) 'During first remapping desc_ac:',lv%desc_ac%is_asb(),&
!!$                 & rmp%desc_ac_pre_remap%is_asb()
!!$            write(0,*) ' First Doing remapping ',lnp, lnp/2
            call psb_remap(lnp/2,rmp%desc_ac_pre_remap,rmp%ac_pre_remap,&
                 & rmp%idest,rmp%isrc,rmp%nrsrc,rmp%naggr,lv%desc_ac,lv%ac,info)
!!$        write(0,*) me,' Out of remapping ',rmp%desc_ac_pre_remap%get_fmt(),' ',&
!!$             & lv%desc_ac%get_fmt(),sum(lv%linmap%naggr),sum(rmp%naggr)
!!$            write(0,*) 'First Assignment ',size(lv%linmap%naggr),size(rmp%naggr)
            lv%linmap%naggr(:) =  rmp%naggr(:)
            lv%linmap%p_desc_V => rmp%desc_ac_pre_remap
            lv%base_a          => lv%ac
            lv%base_desc       => lv%desc_ac
            block
              integer(psb_ipk_) :: meu,npu,mev,npv
              type(psb_ctxt_type) :: ct
              ct = lv%linmap%p_desc_U%get_ctxt()
              call psb_info(ct,meu,npu)
              ct = lv%linmap%p_desc_V%get_ctxt()
              call psb_info(ct,mev,npv)
!!$              write(0,*) 'First Check on out remapping ',i,&
!!$                   & rmp%desc_ac_pre_remap%is_asb(),&
!!$                   & ':',meu,npu,mev,npv
            end block
          end associate
        end if
!!$        write(0,*) 'Second Check on out remapping ',level,&
!!$             & prec%precv(level)%remap_data%desc_ac_pre_remap%is_asb(), newsz
      end if
    end block

    if (info /= psb_success_) then 
      write(ch_err,'(a,i7)') 'Mat asb fail @ level ',i
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err=ch_err)
      goto 9999
    endif
    if (stop_hierarchy_loop) then
      exit array_build_loop
    else
      if (i<iszv) call prec%precv(i)%update_aggr(prec%precv(i+1),info)
    end if
  end do array_build_loop

!!$  write(0,*) ' Done  array_build_loop',iszv,newsz,info,psb_errstatus_fatal()

  if (newsz>0) then
!!$    do i=2,newsz
!!$      write(0,*) me,'Newsz Out of array_build_loop ',i,':',&
!!$           & prec%precv(i)%remap_data%desc_ac_pre_remap%is_asb()
!!$    end do
!!$    write(0,*) 'Calling set_nlevs ',newsz
    call prec%set_nlevs(newsz)
  else
!!$    do i=2, iszv
!!$      write(0,*) me,'Out of array_build_loop ',i,':',&
!!$           & prec%precv(i)%remap_data%desc_ac_pre_remap%is_asb()
!!$    end do
  end if
  iszv = prec%get_nlevs()
  call psb_barrier(ctxt)


!!$  write(0,*) ' Done  reallocating precv',iszv,newsz,info
!!$
!!$  do i=2, iszv
!!$    write(0,*) me,'At end of hierarchy_bld level',i,':',&
!!$         & prec%precv(i)%remap_data%desc_ac_pre_remap%is_asb()
!!$  end do
  call psb_barrier(ctxt)


  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Internal hierarchy build' )
    goto 9999
  endif

  iszv = prec%get_nlevs()
!!$  write(0,*) 'Going for cmp_complexity ',&
!!$       & allocated(prec%precv),iszv,size(prec%precv)
  call prec%cmp_complexity()
  call prec%cmp_avg_cr()

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Exiting with',iszv,' levels'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

#if ( __GNUC__ == 13 && __GNUC_MINOR__ == 3) 
  ! gfortran 13.3.0 generates a strange error here with MOLD
  ! moving to SOURCE but only for this version, since it's heavier
  subroutine save_smoothers(level,save1, save2,info)
    type(amg_d_onelev_type), intent(inout) :: level
    class(amg_d_base_smoother_type), allocatable , intent(inout) :: save1, save2
    integer(psb_ipk_), intent(out) :: info

    info  = 0 
    if (allocated(save1)) then
      call save1%free(info)
      if (info  == 0) deallocate(save1,stat=info)
      if (info /= 0) return
    end if
    if (allocated(save2)) then
      call save2%free(info)
      if (info  == 0) deallocate(save2,stat=info)
      if (info /= 0) return
    end if
    allocate(save1, source=level%sm,stat=info)
    if (info == 0) call level%sm%clone_settings(save1,info)
    if ((info == 0).and.allocated(level%sm2a)) then
      allocate(save2, source=level%sm2a,stat=info) 
      if (info == 0) call level%sm2a%clone_settings(save2,info)
    end if

    return
  end subroutine save_smoothers

  subroutine restore_smoothers(level,save1, save2,info)
    type(amg_d_onelev_type), intent(inout), target :: level
    class(amg_d_base_smoother_type), allocatable, intent(inout) :: save1, save2
    integer(psb_ipk_), intent(out) :: info

    info  = 0

    if (allocated(level%sm)) then
      if (info  == 0) call level%sm%free(info)
      if (info  == 0) deallocate(level%sm,stat=info)
    end if
    if (allocated(save1)) then
      if (info  == 0) allocate(level%sm,source=save1,stat=info)
      if (info == 0) call save1%clone_settings(level%sm,info)
    end if

    if (info /= 0) return

    if (allocated(level%sm2a)) then
      if (info  == 0) call level%sm2a%free(info)
      if (info  == 0) deallocate(level%sm2a,stat=info)
    end if
    if (allocated(save2)) then
      if (info  == 0) allocate(level%sm2a,source=save2,stat=info)
      if (info == 0) call save2%clone_settings(level%sm2a,info)
      if (info  == 0) level%sm2 => level%sm2a
    else
      if (allocated(level%sm)) level%sm2 => level%sm
    end if

    return
  end subroutine restore_smoothers

#else

  subroutine save_smoothers(level,save1, save2,info)
    type(amg_d_onelev_type), intent(inout) :: level
    class(amg_d_base_smoother_type), allocatable , intent(inout) :: save1, save2
    integer(psb_ipk_), intent(out) :: info

    info  = 0 
    if (allocated(save1)) then
      call save1%free(info)
      if (info  == 0) deallocate(save1,stat=info)
      if (info /= 0) return
    end if
    if (allocated(save2)) then
      call save2%free(info)
      if (info  == 0) deallocate(save2,stat=info)
      if (info /= 0) return
    end if
    allocate(save1, mold=level%sm,stat=info)
    if (info == 0) call level%sm%clone_settings(save1,info)
    if ((info == 0).and.allocated(level%sm2a)) then
      allocate(save2, mold=level%sm2a,stat=info) 
      if (info == 0) call level%sm2a%clone_settings(save2,info)
    end if

    return
  end subroutine save_smoothers

  subroutine restore_smoothers(level,save1, save2,info)
    type(amg_d_onelev_type), intent(inout), target :: level
    class(amg_d_base_smoother_type), allocatable, intent(inout) :: save1, save2
    integer(psb_ipk_), intent(out) :: info

    info  = 0

    if (allocated(level%sm)) then
      if (info  == 0) call level%sm%free(info)
      if (info  == 0) deallocate(level%sm,stat=info)
    end if
    if (allocated(save1)) then
      if (info  == 0) allocate(level%sm,mold=save1,stat=info)
      if (info == 0) call save1%clone_settings(level%sm,info)
    end if

    if (info /= 0) return

    if (allocated(level%sm2a)) then
      if (info  == 0) call level%sm2a%free(info)
      if (info  == 0) deallocate(level%sm2a,stat=info)
    end if
    if (allocated(save2)) then
      if (info  == 0) allocate(level%sm2a,mold=save2,stat=info)
      if (info == 0) call save2%clone_settings(level%sm2a,info)
      if (info  == 0) level%sm2 => level%sm2a
    else
      if (allocated(level%sm)) level%sm2 => level%sm
    end if

    return
  end subroutine restore_smoothers
#endif

  function amg_d_policy_do_remap(ctxt,level,aggsize) result(res)
    logical :: res
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: level
    integer(psb_lpk_) :: aggsize
    res =  amg_get_do_remap().and.(level>=2)
!!$    res = .false.
  end function amg_d_policy_do_remap

  subroutine amg_d_hierarchy_bld_cmp_newsz(level,iszv,prevsize,&
       & nlaggr,casize,mnratio,sizeratio,newsz)
    implicit none 
    integer(psb_ipk_) :: level,iszv,newsz
    integer(psb_lpk_) :: nlaggr(:)
    integer(psb_lpk_) :: prevsize, casize
    real(psb_dpk_)    :: mnratio, sizeratio
    ! ==============================
    integer(psb_lpk_) :: iaggsize

    newsz    = 0 
    iaggsize = sum(nlaggr)
    sizeratio = prevsize
    sizeratio = sizeratio/iaggsize
!!$    write(0,*) 'From cmp_newsz: ',iaggsize,casize,&
!!$         & sizeratio,mnratio, level

    if (iaggsize <= casize) newsz = level      
    if (level == iszv)      newsz = level

    if (level>2) then
      if (sizeratio < mnratio) then
        if (sizeratio > 1) then
          newsz = level
        else
          !
          ! We are not gaining 
          !
          newsz = level-1
        end if
      end if
    end if
!!$    write(0,*) 'At end of cmp_newsz ',newsz
  end subroutine amg_d_hierarchy_bld_cmp_newsz

end subroutine amg_d_hierarchy_bld

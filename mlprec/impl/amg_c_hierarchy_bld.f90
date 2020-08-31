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
! File: amg_c_hierarchy_bld.f90
!
! Subroutine: amg_c_hierarchy_bld
! Version:    complex
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
!    a       -  type(psb_cspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_cprec_type), input/output.
!               The preconditioner data structure; upon exit it contains 
!               the multilevel hierarchy of prolongators, restrictors
!               and coarse matrices.
!    info    -  integer, output.
!               Error code.              
!  
subroutine amg_c_hierarchy_bld(a,desc_a,prec,info)

  use psb_base_mod
  use amg_c_inner_mod
  use amg_c_prec_mod, amg_protect_name => amg_c_hierarchy_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  class(amg_cprec_type),intent(inout),target          :: prec
  integer(psb_ipk_), intent(out)                       :: info

  ! Local Variables
  integer(psb_ipk_)  :: ictxt, me,np
  integer(psb_ipk_)  :: err,i,k, err_act, iszv, newsz,&
       & nplevs, mxplevs
  integer(psb_lpk_) :: iaggsize, casize
  real(psb_spk_)     :: mnaggratio, sizeratio, athresh, aomega
  class(amg_c_base_smoother_type), allocatable :: coarse_sm, med_sm, &
       & med_sm2, coarse_sm2
  class(amg_c_base_aggregator_type), allocatable :: tmp_aggr
  type(amg_sml_parms)              :: medparms, coarseparms
  integer(psb_lpk_), allocatable   :: ilaggr(:), nlaggr(:)
  type(psb_lcspmat_type)            :: op_prol
  type(amg_c_onelev_type), allocatable :: tprecv(:)    
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err
  integer(psb_ipk_), save  :: idx_bldtp=-1, idx_matasb=-1
  logical, parameter :: do_timings=.false.

  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'amg_c_hierarchy_bld'
  info = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  prec%ictxt = ictxt
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  if ((do_timings).and.(idx_bldtp==-1))       &
       & idx_bldtp = psb_get_timer_idx("BLD_HIER: bld_tprol")
  if ((do_timings).and.(idx_matasb==-1))     &
       & idx_matasb = psb_get_timer_idx("BLD_HIER: mmat_asb")
  !

  if (.not.allocated(prec%precv)) then 
    !! Error: should have called amg_cprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  newsz      = -1
  mxplevs    = prec%ag_data%max_levs
  mnaggratio = prec%ag_data%min_cr_ratio
  casize     = prec%ag_data%min_coarse_size
  iszv       = size(prec%precv)
  call psb_bcast(ictxt,iszv)
  call psb_bcast(ictxt,casize)
  call psb_bcast(ictxt,mxplevs)
  call psb_bcast(ictxt,mnaggratio)
  if (casize /= prec%ag_data%min_coarse_size) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_coarse_size')
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
  if (iszv /= size(prec%precv)) then 
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
    prec%precv(1)%base_a    => a
    prec%precv(1)%base_desc => desc_a
  
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
  
  if (casize < 0) then
    !
    ! Default to the cubic root of the size at base level.
    ! 
    casize = desc_a%get_global_rows()
    casize = int((sone*casize)**(sone/(sone*3)),psb_lpk_)
    casize = max(casize,lone)
    casize = casize*40_psb_lpk_
    call psb_bcast(ictxt,casize)
    if (casize > huge(prec%ag_data%min_coarse_size)) then
      !
      ! computed coarse size does not fit in IPK_.
      ! This is very unlikely, but make sure to put a positive number
      !
      prec%ag_data%min_coarse_size = huge(prec%ag_data%min_coarse_size)
    else
      prec%ag_data%min_coarse_size = casize
    end if
  end if
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
  ! First set desired number of levels
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
    iszv = size(prec%precv)
  end if

  !
  ! Finest level first; create a GEN_BLOCK
  ! copy of the descriptor. 
  ! 
  prec%precv(1)%base_a    => a
  call psb_cd_renum_block(desc_a,prec%precv(1)%desc_ac,info)
  prec%precv(1)%base_desc => prec%precv(1)%desc_ac
  newsz = 0
  array_build_loop: do i=2, iszv
    !
    ! Check on the iprcparm contents: they should be the same
    ! on all processes.
    !
    call psb_bcast(ictxt,prec%precv(i)%parms)

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
    ! Build the mapping between levels i-1 and i and the matrix
    ! at level i
    !
    if (do_timings) call psb_tic(idx_bldtp)
    if (info == psb_success_)&
         & call prec%precv(i)%bld_tprol(prec%precv(i-1)%base_a,&
         & prec%precv(i-1)%base_desc,&
         & ilaggr,nlaggr,op_prol,prec%ag_data,info)
    if (do_timings) call psb_toc(idx_bldtp)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Map build')
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
    iaggsize = sum(nlaggr)

    sizeratio = iaggsize
    if (i==2) then 
      sizeratio = desc_a%get_global_rows()/sizeratio
    else
      sizeratio = sum(prec%precv(i-1)%map%naggr)/sizeratio
    end if
    prec%precv(i)%szratio = sizeratio

    if (iaggsize <= casize) newsz = i      
    if (i == iszv)          newsz = i

    if (i>2) then
      if (sizeratio < mnaggratio) then
        if (sizeratio > 1) then
          newsz = i
        else
          !
          ! We are not gaining 
          !
          newsz = i-1
        end if
      end if

      if (all(nlaggr == prec%precv(i-1)%map%naggr)) then 
        newsz=i-1
        if (me == 0) then 
          write(debug_unit,*) trim(name),&
               &': Warning: aggregates from level ',&
               & newsz
          write(debug_unit,*) trim(name),&
               &':                       to level ',&
               & iszv,' coincide.'
          write(debug_unit,*) trim(name),&
               &': Number of levels actually used :',newsz
          write(debug_unit,*)
        end if
      end if
    end if
    call psb_bcast(ictxt,newsz)
    
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
        ! We are going back and revisit a previous leve;
        ! recover the aggregation.
        !
        ilaggr = prec%precv(newsz)%map%iaggr
        nlaggr = prec%precv(newsz)%map%naggr
        call prec%precv(newsz)%tprol%clone(op_prol,info)
      end if
      if (do_timings) call psb_tic(idx_matasb)      
      if (info == psb_success_) call prec%precv(newsz)%mat_asb( &
           & prec%precv(newsz-1)%base_a,prec%precv(newsz-1)%base_desc,&
           & ilaggr,nlaggr,op_prol,info)
      if (do_timings) call psb_toc(idx_matasb)      
      if (info /= 0) then 
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Mat asb')
        goto 9999
      endif
      exit array_build_loop
    else
      if (do_timings) call psb_tic(idx_matasb)            
      if (info == psb_success_) call prec%precv(i)%mat_asb(&
           & prec%precv(i-1)%base_a,prec%precv(i-1)%base_desc,&
           & ilaggr,nlaggr,op_prol,info)
      if (do_timings) call psb_toc(idx_matasb)      
    end if
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Map build')
      goto 9999
    endif
    if (i<iszv) call prec%precv(i)%update_aggr(prec%precv(i+1),info)

  end do array_build_loop

  if (newsz > 0) then
    !
    ! We exited early from the build loop, need to fix
    ! the size.
    !
    allocate(tprecv(newsz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=1,newsz
      call prec%precv(i)%move_alloc(tprecv(i),info)
    end do
    do i=newsz+1, iszv
      call prec%precv(i)%free(info)
    end do
    call move_alloc(tprecv,prec%precv) 
    ! Ignore errors from transfer
    info = psb_success_
    !
    ! Restart
    iszv = newsz
    ! Fix the pointers, but the level 1 should
    ! be treated differently
    if (.not.associated(prec%precv(1)%base_desc,desc_a)) then
      prec%precv(1)%base_desc => prec%precv(1)%desc_ac
    end if
    do i=2, iszv 
      prec%precv(i)%base_a       => prec%precv(i)%ac
      prec%precv(i)%base_desc    => prec%precv(i)%desc_ac
      prec%precv(i)%map%p_desc_U => prec%precv(i-1)%base_desc
      prec%precv(i)%map%p_desc_V => prec%precv(i)%base_desc
    end do
  end if

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Internal hierarchy build' )
    goto 9999
  endif

  iszv = size(prec%precv)

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
  subroutine save_smoothers(level,save1, save2,info)
    type(amg_c_onelev_type), intent(inout) :: level
    class(amg_c_base_smoother_type), allocatable , intent(inout) :: save1, save2
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
    type(amg_c_onelev_type), intent(inout), target :: level
    class(amg_c_base_smoother_type), allocatable, intent(inout) :: save1, save2
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
  
end subroutine amg_c_hierarchy_bld

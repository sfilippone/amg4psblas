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
! File: amg_c_extprol_bld.f90
!
! Subroutine: amg_c_extprol_bld
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
!    a       -  type(psb_cspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_cprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!
!    amold   -  class(psb_c_base_sparse_mat), input, optional
!               Mold for the inner format of matrices contained in the
!               preconditioner
!
!
!    vmold   -  class(psb_c_base_vect_type), input, optional
!               Mold for the inner format of vectors contained in the
!               preconditioner
!
!
!  
subroutine amg_c_extprol_bld(a,desc_a,p,prolv,restrv,info,amold,vmold,imold)

  use psb_base_mod
  use amg_c_inner_mod
  use amg_c_prec_mod, amg_protect_name => amg_c_extprol_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_cspmat_type),intent(inout), target        :: prolv(:)
  type(psb_cspmat_type),intent(inout), target        :: restrv(:)
  type(psb_desc_type), intent(inout), target         :: desc_a
  type(amg_cprec_type),intent(inout),target          :: p
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_c_base_sparse_mat), intent(in), optional :: amold
  class(psb_c_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  ! !$  character, intent(in), optional         :: upd

  ! Local Variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me,np
  integer(psb_ipk_)   :: err,i,k, err_act, iszv, newsz, casize, nplevs, mxplevs
  integer(psb_ipk_)   :: nprolv, nrestrv
  real(psb_spk_)      :: mnaggratio
  integer(psb_ipk_)  :: ipv(amg_ifpsz_), val
  class(amg_c_base_smoother_type), allocatable :: coarse_sm, base_sm, med_sm
  type(amg_sml_parms)              :: baseparms, medparms, coarseparms
  type(amg_c_onelev_type), allocatable :: tprecv(:)    
  integer(psb_ipk_)  :: int_err(5)
  character          :: upd_
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err
  logical, parameter :: debug=.false.

  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'amg_c_extprol_bld'
  info = psb_success_
  int_err(1) = 0
  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  p%ctxt = ctxt
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
#if defined(LPK8)
  info=psb_err_internal_error_
  call psb_errpush(info,name,a_err='Need fix for LPK8')
  goto 9999
#else
  
  !
  ! For the time being we are commenting out the UPDATE argument
  ! we plan to resurrect it later. 
  ! !$  if (present(upd)) then 
  ! !$    if (debug_level >= psb_debug_outer_) &
  ! !$         & write(debug_unit,*) me,' ',trim(name),'UPD ', upd
  ! !$
  ! !$    if ((psb_toupper(upd).eq.'F').or.(psb_toupper(upd).eq.'T')) then
  ! !$      upd_=psb_toupper(upd)
  ! !$    else
  ! !$      upd_='F'
  ! !$    endif
  ! !$  else
  ! !$    upd_='F'
  ! !$  endif
  upd_ = 'F'

  if (.not.allocated(p%precv)) then 
    !! Error: should have called amg_cprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if
  
  !
  ! Check to ensure all procs have the same 
  !   
  newsz      = -1
  mxplevs    = p%ag_data%max_levs
  mnaggratio = p%ag_data%min_cr_ratio
  casize     = p%ag_data%min_coarse_size
  iszv       = size(p%precv)
  nprolv     = size(prolv)
  nrestrv    = size(restrv)
  call psb_bcast(ctxt,iszv)
  call psb_bcast(ctxt,casize)
  call psb_bcast(ctxt,mxplevs)
  call psb_bcast(ctxt,mnaggratio)
  call psb_bcast(ctxt,nprolv)
  call psb_bcast(ctxt,nrestrv)
  if (casize /= p%ag_data%min_coarse_size) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_coarse_size')
    goto 9999
  end if
  if (mxplevs /= p%ag_data%max_levs) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent max_levs')
    goto 9999
  end if
  if (mnaggratio /= p%ag_data%min_cr_ratio) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent min_cr_ratio')
    goto 9999
  end if
  if (iszv /= size(p%precv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if
  if (nprolv /= size(prolv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of prolv')
    goto 9999
  end if
  if (nrestrv /= size(restrv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of restrv')
    goto 9999
  end if
  if (nrestrv /= nprolv) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size prolv vs restrv')
    goto 9999
  end if

  if (iszv <= 1) then
    ! We should only ever get here for multilevel.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  if (nrestrv < 1) then
    ! We should only ever get here for multilevel.
    info=psb_err_from_subroutine_
    ch_err='size restrv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  !
  nplevs     =  nrestrv + 1
  p%ag_data%max_levs = nplevs

  ! 
  !  Fixed number of levels. 
  !
  nplevs = max(itwo,mxplevs)

  coarseparms = p%precv(iszv)%parms
  baseparms   = p%precv(1)%parms
  medparms    = p%precv(2)%parms

  allocate(coarse_sm, source=p%precv(iszv)%sm,stat=info) 
  if (info == psb_success_)  &
       & allocate(med_sm, source=p%precv(2)%sm,stat=info) 
  if (info == psb_success_)  &
       & allocate(base_sm, source=p%precv(1)%sm,stat=info) 
  if (info /= psb_success_) then 
    write(0,*) 'Error in saving smoothers',info
    call psb_errpush(psb_err_internal_error_,name,a_err='Base level precbuild.')
    goto 9999
  end if

  if (iszv /= nplevs) then
    allocate(tprecv(nplevs),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    tprecv(1)%parms = baseparms
    allocate(tprecv(1)%sm,source=base_sm,stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=2,nplevs-1
      tprecv(i)%parms = medparms
      allocate(tprecv(i)%sm,source=med_sm,stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,&
             & a_err='prec reallocation')
        goto 9999
      endif
    end do
    tprecv(nplevs)%parms = coarseparms
    allocate(tprecv(nplevs)%sm,source=coarse_sm,stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,&
           & a_err='prec reallocation')
      goto 9999
    endif
    do i=1,iszv
      call p%precv(i)%free(info)
    end do
    call move_alloc(tprecv,p%precv)
    iszv = size(p%precv)
  end if
  !
  ! Finest level first; remember to fix base_a and base_desc
  ! 
  p%precv(1)%base_a    => a
  p%precv(1)%base_desc => desc_a
  newsz = 0 
  array_build_loop: do i=2, iszv

    !
    ! Sanity checks on the parameters
    !
    if (i<iszv) then 
      !
      ! A replicated matrix only makes sense at the coarsest level
      !
      call amg_check_def(p%precv(i)%parms%coarse_mat,'Coarse matrix',&
           &   amg_distr_mat_,is_distr_ml_coarse_mat)
    end if
    if (debug.and.(me==0)) write(0,*)name,' Building aggregation at level ',i
    call amg_c_extaggr_bld(p%precv(i-1)%base_a,&
         & p%precv(i-1)%base_desc,p%precv(i),restrv(i-1),prolv(i-1),info)
    p%precv(i)%base_a    => p%precv(i)%ac
    p%precv(i)%base_desc => p%precv(i)%desc_ac

    
    if (i>2) then 
      if (all(p%precv(i)%map%naggr == p%precv(i-1)%map%naggr)) then 
        newsz=i-1
      end if
      call psb_bcast(ctxt,newsz)
      if (newsz > 0) exit array_build_loop
    end if
  end do array_build_loop

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Internal extprol build' )
    goto 9999
  endif

  iszv = size(p%precv)
  
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Exiting with',iszv,' levels'
#endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
contains

  subroutine amg_c_extaggr_bld(a,desc_a,p,op_restr,op_prol,info)
    use psb_base_mod
    use amg_c_inner_mod    

    implicit none

    ! Arguments
    type(psb_cspmat_type), intent(in), target     :: a
    type(psb_cspmat_type), intent(inout)     :: op_restr,op_prol
    type(psb_desc_type), intent(in), target       :: desc_a
    type(amg_c_onelev_type), intent(inout),target :: p
    integer(psb_ipk_), intent(out)                :: info

    ! Local variables
    character(len=20)   :: name
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_)   :: np, me, ncol
    integer(psb_ipk_)   :: err_act,ntaggr,nzl 
    integer(psb_ipk_), allocatable   :: ilaggr(:), nlaggr(:)
    type(psb_cspmat_type)      :: ac, am2, am3, am4
    type(psb_c_coo_sparse_mat) :: acoo, bcoo
    type(psb_c_csr_sparse_mat) :: acsr1
    logical, parameter :: debug=.false.

    name='amg_c_extaggr_bld'
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if
    info = psb_success_
    ctxt = desc_a%get_context()
    call psb_info(ctxt,me,np)
#if defined(LPK8)
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Need fix for LPK8')
    goto 9999
#else
    allocate(nlaggr(np),ilaggr(1))
    nlaggr = 0
    ilaggr = 0
    p%parms%par_aggr_alg = amg_ext_aggr_
    call amg_check_def(p%parms%ml_cycle,'Multilevel cycle',&
         &   amg_mult_ml_,is_legal_ml_cycle)
    call amg_check_def(p%parms%coarse_mat,'Coarse matrix',&
         &   amg_distr_mat_,is_legal_ml_coarse_mat)

    nlaggr(me+1) = op_restr%get_nrows()
    if (op_restr%get_nrows() /= op_prol%get_ncols()) then
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='Inconsistent restr/prol sizes')
      goto 9999      
    end if
    call psb_sum(ctxt,nlaggr)
    ntaggr = sum(nlaggr)
    ncol = desc_a%get_local_cols()
    if (debug) write(0,*)me,' Sizes:',op_restr%get_nrows(),op_restr%get_ncols(),&
         & op_prol%get_nrows(),op_prol%get_ncols(), a%get_nrows(),a%get_ncols()
    !
    ! Compute local part of AC
    !
    call op_prol%clone(am2,info)
    if (info == psb_success_) call psb_sphalo(am2,desc_a,am4,info,&
         & colcnv=.false.,rowscale=.true.)
    if (info == psb_success_) call psb_rwextd(ncol,am2,info,b=am4)
    if (info == psb_success_) call am4%free()
    call psb_spspmm(a,am2,am3,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 2')
      goto 9999
    end if
    call psb_sphalo(am3,desc_a,am4,info,&
         & colcnv=.false.,rowscale=.true.)
    if (info == psb_success_) call psb_rwextd(ncol,am3,info,b=am4)      
    if (info == psb_success_) call am4%free()
    if(info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
      goto 9999
    end if
    call psb_spspmm(op_restr,am3,ac,info)
    if (info == psb_success_) call am3%free()
    if (info == psb_success_) call ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,a_err='Build ac = op_restr x am3')
      goto 9999
    end if
    
    select case(p%parms%coarse_mat)

    case(amg_distr_mat_) 

      call ac%mv_to(bcoo)
      nzl = bcoo%get_nzeros()

      if (info == psb_success_) call psb_cdall(ctxt,p%desc_ac,info,nl=nlaggr(me+1))
      if (info == psb_success_) call psb_cdins(nzl,bcoo%ia,bcoo%ja,p%desc_ac,info)
      if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
      if (info == psb_success_) call psb_glob_to_loc(bcoo%ia(1:nzl),p%desc_ac,info,iact='I')
      if (info == psb_success_) call psb_glob_to_loc(bcoo%ja(1:nzl),p%desc_ac,info,iact='I')
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Creating p%desc_ac and converting ac')
        goto 9999
      end if
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Assembld aux descr. distr.'
      call p%ac%mv_from(bcoo)

      call p%ac%set_nrows(p%desc_ac%get_local_rows())
      call p%ac%set_ncols(p%desc_ac%get_local_cols())
      call p%ac%set_asb()

      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
        goto 9999
      end if

      if (np>1) then 
        call op_prol%mv_to(acsr1)
        nzl = acsr1%get_nzeros()
        call psb_glob_to_loc(acsr1%ja(1:nzl),p%desc_ac,info,'I')
        if(info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_glob_to_loc')
          goto 9999
        end if
        call op_prol%mv_from(acsr1)
      endif
      call op_prol%set_ncols(p%desc_ac%get_local_cols())

      if (np>1) then 
        call op_restr%cscnv(info,type='coo',dupl=psb_dupl_add_)
        call op_restr%mv_to(acoo)
        nzl = acoo%get_nzeros()
        if (info == psb_success_) call psb_glob_to_loc(acoo%ia(1:nzl),p%desc_ac,info,'I')
        call acoo%set_dupl(psb_dupl_add_)
        if (info == psb_success_) call op_restr%mv_from(acoo)
        if (info == psb_success_) call op_restr%cscnv(info,type='csr')        
        if(info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Converting op_restr to local')
          goto 9999
        end if
      end if
      call op_restr%set_nrows(p%desc_ac%get_local_cols())

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & 'Done ac '

    case(amg_repl_mat_) 
      !
      !
      call psb_cdall(ctxt,p%desc_ac,info,mg=ntaggr,repl=.true.)
      if (info == psb_success_) call psb_cdasb(p%desc_ac,info)
      if (info == psb_success_) &
           & call psb_gather(p%ac,ac,p%desc_ac,info,dupl=psb_dupl_add_,keeploc=.false.)

      if (info /= psb_success_) goto 9999

    case default 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid amg_coarse_mat_')
      goto 9999
    end select

    call p%ac%cscnv(info,type='csr',dupl=psb_dupl_add_)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
      goto 9999
    end if
    
    !
    ! Copy the prolongation/restriction matrices into the descriptor map.
    !  op_restr => PR^T   i.e. restriction  operator
    !  op_prol => PR     i.e. prolongation operator
    !  
    
    p%map = psb_linmap(psb_map_aggr_,desc_a,&
         & p%desc_ac,op_restr,op_prol,ilaggr,nlaggr)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_Free')
      goto 9999
    end if
#endif
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine amg_c_extaggr_bld
    
end subroutine amg_c_extprol_bld

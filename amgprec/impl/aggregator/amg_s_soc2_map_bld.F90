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
!
! File: amg_s_soc2_map__bld.f90
!
! Subroutine: amg_s_soc2_map_bld
! Version:    real
!
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. This variant is based on the method
!  presented in 
!
!    S. Gratton, P. Henon, P. Jiranek and X. Vasseur:
!    Reducing complexity of algebraic multigrid by aggregation
!    Numerical Lin. Algebra with Applications, 2016, 23:501-518
!
! Note: upon exit 
!
! Arguments:
!    a       -  type(psb_sspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_sprec_type), input/output.
!               The preconditioner data structure; upon exit it contains 
!               the multilevel hierarchy of prolongators, restrictors
!               and coarse matrices.
!    info    -  integer, output.
!               Error code.              
!
!
!
subroutine amg_s_soc2_map_bld(iorder,theta,clean_zeros,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod 
  use amg_base_prec_type
  use amg_s_inner_mod
#if defined(OPENMP)
  use omp_lib
#endif

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)     :: iorder
  logical, intent(in)               :: clean_zeros
  type(psb_sspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)    :: desc_a
  real(psb_spk_), intent(in)         :: theta
  integer(psb_lpk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
  integer(psb_ipk_), intent(out)               :: info

  ! Local variables
  integer(psb_ipk_), allocatable  :: ils(:), neigh(:), irow(:), icol(:),&
       & ideg(:), idxs(:)
  integer(psb_lpk_), allocatable :: tmpaggr(:)
  real(psb_spk_), allocatable  :: val(:), diag(:)
  integer(psb_ipk_) :: icnt,nlp,k,n,ia,isz,nr,nc,naggr,i,j,m, nz, ilg, ii, ip, ip1,nzcnt
  integer(psb_lpk_) :: nrglob
  type(psb_s_csr_sparse_mat) :: acsr, muij, s_neigh
  type(psb_s_coo_sparse_mat) :: s_neigh_coo
  real(psb_spk_)  :: cpling, tcl
  logical :: disjoint
  integer(psb_ipk_) :: debug_level, debug_unit,err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me
  integer(psb_ipk_) :: nrow, ncol, n_ne
  character(len=20)  :: name, ch_err
  integer(psb_ipk_), save :: idx_soc2_p1=-1, idx_soc2_p2=-1, idx_soc2_p3=-1
  integer(psb_ipk_), save :: idx_soc2_p0=-1
  logical, parameter      :: do_timings=.true.

  info=psb_success_
  name = 'amg_soc2_map_bld'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ctxt=desc_a%get_context()
  call psb_info(ctxt,me,np)
  nrow   = desc_a%get_local_rows()
  ncol   = desc_a%get_local_cols()
  nrglob = desc_a%get_global_rows()
  if ((do_timings).and.(idx_soc2_p0==-1))       &
       & idx_soc2_p0 = psb_get_timer_idx("SOC2_MAP: phase0")
  if ((do_timings).and.(idx_soc2_p1==-1))       &
       & idx_soc2_p1 = psb_get_timer_idx("SOC2_MAP: phase1")
  if ((do_timings).and.(idx_soc2_p2==-1))       &
       & idx_soc2_p2 = psb_get_timer_idx("SOC2_MAP: phase2")
  if ((do_timings).and.(idx_soc2_p3==-1))       &
       & idx_soc2_p3 = psb_get_timer_idx("SOC2_MAP: phase3")

  nr = a%get_nrows()
  nc = a%get_ncols()
  allocate(ilaggr(nr),neigh(nr),ideg(nr),idxs(nr),icol(nc),stat=info)
  if(info /= psb_success_) then
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/2*nr,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if

  if (do_timings) call psb_tic(idx_soc2_p0)
  diag = a%get_diag(info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_sp_getdiag')
    goto 9999
  end if

  !
  ! Phase zero: compute muij
  ! 
  call a%cp_to(muij)
  if (clean_zeros) call muij%clean_zeros(info)
  !$omp parallel do private(i,j,k) shared(nr,diag,muij) schedule(static)
  do i=1, nr
    do k=muij%irp(i),muij%irp(i+1)-1
      j = muij%ja(k)
      if (j<= nr) muij%val(k) = abs(muij%val(k))/sqrt(abs(diag(i)*diag(j)))
    end do
  end do
  !$omp end parallel do 
  !
  ! Compute the 1-neigbour; mark strong links with +1, weak links with -1
  !
  call s_neigh_coo%allocate(nr,nr,muij%get_nzeros())
  !$omp parallel do private(i,j,k) shared(nr,diag,muij) schedule(static)
  do i=1, nr
    do k=muij%irp(i),muij%irp(i+1)-1
      j = muij%ja(k)
      s_neigh_coo%ia(k)  = i
      s_neigh_coo%ja(k)  = j
      if (j<=nr) then 
        if (real(muij%val(k)) >= theta) then 
          s_neigh_coo%val(k) = sone
        else
          s_neigh_coo%val(k) = -sone
        end if
      else
        s_neigh_coo%val(k) = -sone        
      end if
    end do
  end do
  !$omp end parallel do
  !write(*,*) 'S_NEIGH: ',nr,ip
  call s_neigh_coo%set_nzeros(muij%get_nzeros())
  call s_neigh%mv_from_coo(s_neigh_coo,info)

  if (iorder == amg_aggr_ord_nat_) then
    
    !$omp parallel do private(i) shared(ilaggr,idxs) schedule(static)
    do i=1, nr
      ilaggr(i) = -(nr+1)
      idxs(i)   = i 
    end do
    !$omp end parallel do 
  else 
    !$omp parallel do private(i) shared(ilaggr,idxs,muij)  schedule(static)
    do i=1, nr
      ilaggr(i) = -(nr+1)
      ideg(i)   = muij%irp(i+1) - muij%irp(i)
    end do
    !$omp end parallel do
    call psb_msort(ideg,ix=idxs,dir=psb_sort_down_)
  end if

  if (do_timings) call psb_toc(idx_soc2_p0)
  if (do_timings) call psb_tic(idx_soc2_p1)

  !
  ! Phase one: Start with disjoint groups.
  ! 
  naggr = 0
#if defined(OPENMP)
  block
    integer(psb_ipk_), allocatable  :: bnds(:), locnaggr(:)
    integer(psb_ipk_) :: myth,nths, kk
    ! The parallelization makes use of a locaggr(:) array; each thread
    ! keeps its own version of naggr, and when the loop ends, a prefix is applied
    ! to locnaggr to determine:
    ! 1. The total number of aggregaters NAGGR;
    ! 2. How much should each thread shift its own aggregates
    ! Part 2 requires to keep track of which thread defined each entry
    ! of ilaggr(), so that each entry can be adjusted correctly: even
    ! if an entry I belongs to the range BNDS(TH)>BNDS(TH+1)-1, it may have
    ! been set because it is strongly connected to an entry J belonging to a
    ! different thread. 

    !$omp parallel shared(s_neigh,bnds,idxs,locnaggr,ilaggr,nr,naggr,diag,theta,nths,info) &
    !$omp private(icol,val,myth,kk) 
    block
      integer(psb_ipk_) :: ii,nlp,k,kp,n,ia,isz,nc,i,j,m,nz,ilg,ip,rsz,ip1,nzcnt
      integer(psb_lpk_) :: itmp
      !$omp master
      nths = omp_get_num_threads()
      allocate(bnds(0:nths),locnaggr(0:nths+1))
      locnaggr(:) = 0
      bnds(0) = 1
      !$omp end master      
      !$omp barrier
      myth = omp_get_thread_num()
      rsz = nr/nths
      if (myth < mod(nr,nths)) rsz = rsz + 1
      bnds(myth+1) = rsz
      !$omp barrier
      !$omp master
      do i=1,nths
        bnds(i) = bnds(i) + bnds(i-1)
      end do
      info = 0
      !$omp end master
      !$omp barrier

      !$omp do schedule(static) private(disjoint) 
      do kk=0, nths-1
        step1: do ii=bnds(kk), bnds(kk+1)-1
          i = idxs(ii)
          if (info /= 0) then
            write(0,*) ' Step1:',kk,ii,i,info
            cycle step1
          end if
          if ((i<1).or.(i>nr)) then
            !$omp atomic write
            info=psb_err_internal_error_
            !$omp end atomic 
            call psb_errpush(info,name)
            cycle step1
            !goto 9999
          end if


          if (ilaggr(i) == -(nr+1)) then 
            !
            ! Get the 1-neighbourhood of I 
            !
            ip1 = s_neigh%irp(i)
            nz  = s_neigh%irp(i+1)-ip1
            !
            ! If the neighbourhood only contains I, skip it
            !
            if (nz ==0) then
              ilaggr(i) = 0
              cycle step1
            end if
            if ((nz==1).and.(s_neigh%ja(ip1)==i)) then
              ilaggr(i) = 0
              cycle step1
            end if

            nzcnt = count(real(s_neigh%val(ip1:ip1+nz-1)) > 0)
            icol(1:nzcnt) = pack(s_neigh%ja(ip1:ip1+nz-1),(real(s_neigh%val(ip1:ip1+nz-1)) > 0))
            disjoint = all(ilaggr(icol(1:nzcnt)) == -(nr+1)) 

            !
            ! If the whole strongly coupled neighborhood of I is
            ! as yet unconnected, turn it into the next aggregate.
            ! Same if ip==0 (in which case, neighborhood only
            ! contains I even if it does not look like it from matrix)
            ! The fact that DISJOINT is private and not under lock
            ! generates a certain un-repeatability, in that between
            ! computing DISJOINT and assigning, another thread might
            ! alter the values of ILAGGR.
            ! However, a certain unrepeatability is already present
            ! because the sequence of aggregates is computed with a
            ! different order than in serial mode.
            ! In any case, even if the enteries of ILAGGR may be
            ! overwritten, the important thing is that each entry is
            ! consistent and they generate a correct aggregation map.
            !
            if (disjoint) then
              locnaggr(kk)     = locnaggr(kk) + 1
              itmp = (bnds(kk)-1+locnaggr(kk))*nths+kk
              if (itmp < (bnds(kk)-1+locnaggr(kk))) then
                !$omp atomic update
                info = max(12345678,info)
                !$omp end atomic
                cycle step1
              end if
              !$omp atomic write
              ilaggr(i) = itmp
              !$omp end atomic
              do k=1, nzcnt
                !$omp atomic write
                ilaggr(icol(k)) = itmp
                !$omp end atomic
              end do
            end if
          end if
        enddo step1
      end do
      !$omp end do

      !$omp master
      naggr = sum(locnaggr(0:nths-1))
      do i=1,nths
        locnaggr(i) = locnaggr(i) + locnaggr(i-1)
      end do
      do i=nths+1,1,-1
        locnaggr(i) = locnaggr(i-1)
      end do
      locnaggr(0) = 0
      !write(0,*) 'LNAG ',locnaggr(nths+1)
      !$omp end master 
      !$omp barrier
      !$omp  do schedule(static) 
      do kk=0, nths-1
        do ii=bnds(kk), bnds(kk+1)-1
          if (ilaggr(ii) > 0) then 
            kp = mod(ilaggr(ii),nths)
            ilaggr(ii) = (ilaggr(ii)/nths)- (bnds(kp)-1) + locnaggr(kp)
          end if
        end do
      end do
      !$omp end do
    end block
    !$omp end parallel
  end block
  if (info /= 0) then
    if (info == 12345678) write(0,*) 'Overflow in encoding ILAGGR'
    info=psb_err_internal_error_
    call psb_errpush(info,name)
    goto 9999
  end if

#else
  icnt = 0
  step1: do ii=1, nr
    i = idxs(ii)

    if (ilaggr(i) == -(nr+1)) then 
      !
      ! Get the 1-neighbourhood of I 
      !
      ip1 = s_neigh%irp(i)
      nz  = s_neigh%irp(i+1)-ip1
      !
      ! If the neighbourhood only contains I, skip it
      !
      if (nz ==0) then
        ilaggr(i) = 0
        cycle step1
      end if
      if ((nz==1).and.(s_neigh%ja(ip1)==i)) then
        ilaggr(i) = 0
        cycle step1
      end if      
      !
      ! If the whole strongly coupled neighborhood of I is
      ! as yet unconnected, turn it into the next aggregate.
      !
      nzcnt = count(real(s_neigh%val(ip1:ip1+nz-1)) > 0)
      icol(1:nzcnt) = pack(s_neigh%ja(ip1:ip1+nz-1),(real(s_neigh%val(ip1:ip1+nz-1)) > 0))
      disjoint = all(ilaggr(icol(1:nzcnt)) == -(nr+1)) 
      if (disjoint) then 
        icnt      = icnt + 1 
        naggr     = naggr + 1
        do k=1, nzcnt
          ilaggr(icol(k)) = naggr
        end do
        ilaggr(i) = naggr
      end if
    endif
  enddo step1
#endif  
  if (debug_level >= psb_debug_outer_) then 
    write(debug_unit,*) me,' ',trim(name),&
         & ' Check 1:',count(ilaggr == -(nr+1))
  end if

  !
  ! Phase two: join the neighbours
  !
  tmpaggr = ilaggr
  step2: do ii=1,nr
    i = idxs(ii)

    if (ilaggr(i) == -(nr+1)) then         
      !
      ! Find the most strongly connected neighbour that is
      ! already aggregated, if any, and join its aggregate
      !
      cpling = szero
      ip = 0
      do k=s_neigh%irp(i), s_neigh%irp(i+1)-1
        j   = s_neigh%ja(k)
        if ((1<=j).and.(j<=nr)) then 
          if ( (tmpaggr(j) > 0).and. (real(muij%val(k)) > cpling)&
               & .and.(real(s_neigh%val(k))>0)) then
            ip = k
            cpling = muij%val(k)
          end if
        end if
      enddo
      if (ip > 0) then 
        ilaggr(i) = ilaggr(s_neigh%ja(ip))
      end if
    end if
  end do step2


  !
  ! Phase three: sweep over leftovers, if any 
  !
  step3: do ii=1,nr
    i = idxs(ii)

    if (ilaggr(i) < 0) then
      !
      ! Find its strongly  connected neighbourhood not 
      ! already aggregated, and make it into a new aggregate.
      !
      ip = 0 
      do k=s_neigh%irp(i), s_neigh%irp(i+1)-1
        j   = s_neigh%ja(k)
        if ((1<=j).and.(j<=nr)) then 
          if (ilaggr(j) < 0)  then
            ip       = ip + 1
            icol(ip) = j
          end if
        end if
      enddo
      if (ip > 0) then
        icnt      = icnt + 1 
        naggr     = naggr + 1
        ilaggr(i) = naggr
        do k=1, ip
          ilaggr(icol(k)) = naggr
        end do
      end if
    end if
  end do step3

  ! Any leftovers?
  do i=1, nr
    if (ilaggr(i) <= 0) then
      nz = (s_neigh%irp(i+1)-s_neigh%irp(i))
      if (nz <= 1) then
        ! Mark explicitly as a singleton so that 
        ! it will be ignored in map_to_tprol.
        ! Need to use -(nrglob+nr) to make sure
        ! it's still negative when shifted and combined with
        ! other processes. 
        ilaggr(i) = -(nrglob+nr)
      else
        info=psb_err_internal_error_
        call psb_errpush(info,name,a_err='Fatal error: non-singleton leftovers')
        goto 9999
      endif
    end if
  end do

  if (naggr > ncol) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Fatal error: naggr>ncol')
    goto 9999
  end if

  call psb_realloc(ncol,ilaggr,info)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  allocate(nlaggr(np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/np,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if

  nlaggr(:) = 0
  nlaggr(me+1) = naggr
  call psb_sum(ctxt,nlaggr(1:np))

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_s_soc2_map_bld


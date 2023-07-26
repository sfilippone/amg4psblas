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
! File: amg_d_soc1_map__bld.f90
!
! Subroutine: amg_d_soc1_map_bld
! Version:    real
!
!  This routine builds the tentative prolongator based on the
!  strength of connection aggregation algorithm presented in
!
!    M. Brezina and P. Vanek, A black-box iterative solver based on a 
!    two-level Schwarz method, Computing,  63 (1999), 233-263.
!    P. Vanek, J. Mandel and M. Brezina, Algebraic Multigrid by Smoothed
!    Aggregation for Second and Fourth Order Elliptic Problems, Computing, 56
!    (1996), 179-196.
!
! Note: upon exit 
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
!
!
subroutine amg_d_soc1_map_bld(iorder,theta,clean_zeros,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod
  use amg_base_prec_type
  use amg_d_inner_mod
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)     :: iorder
  logical, intent(in)               :: clean_zeros
  type(psb_dspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)    :: desc_a
  real(psb_dpk_), intent(in)         :: theta
  integer(psb_lpk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
  integer(psb_ipk_), intent(out)               :: info

  ! Local variables
  integer(psb_ipk_), allocatable :: ioffs(:), neigh(:), irow(:), icol(:),&
       & ideg(:), idxs(:)
  integer(psb_lpk_), allocatable :: tmpaggr(:)
  real(psb_dpk_), allocatable  :: val(:), diag(:)
  integer(psb_ipk_) :: icnt,nlp,k,n,ia,isz,nr, nc, naggr,i,j,m, nz, ilg, ii, ip
  type(psb_d_csr_sparse_mat) :: acsr
  real(psb_dpk_)  :: cpling, tcl
  logical :: disjoint
  integer(psb_ipk_) :: debug_level, debug_unit,err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_ipk_)   :: nrow, ncol, n_ne
  integer(psb_lpk_)   :: nrglob
  character(len=20)   :: name, ch_err
  integer(psb_ipk_), save :: idx_soc1_p1=-1, idx_soc1_p2=-1, idx_soc1_p3=-1
  integer(psb_ipk_), save :: idx_soc1_p0=-1
  logical, parameter      :: do_timings=.true.

  info=psb_success_
  name = 'amg_soc1_map_bld'
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
  if ((do_timings).and.(idx_soc1_p0==-1))       &
       & idx_soc1_p0 = psb_get_timer_idx("SOC1_MAP: phase0")
  if ((do_timings).and.(idx_soc1_p1==-1))       &
       & idx_soc1_p1 = psb_get_timer_idx("SOC1_MAP: phase1")
  if ((do_timings).and.(idx_soc1_p2==-1))       &
       & idx_soc1_p2 = psb_get_timer_idx("SOC1_MAP: phase2")
  if ((do_timings).and.(idx_soc1_p3==-1))       &
       & idx_soc1_p3 = psb_get_timer_idx("SOC1_MAP: phase3")

  nr = a%get_nrows()
  nc = a%get_ncols()
  allocate(ilaggr(nr),ioffs(nr),neigh(nr),ideg(nr),idxs(nr),&
       & icol(nc),val(nc),stat=info)
  if(info /= psb_success_) then
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/2*nr,izero,izero,izero,izero/),&
         & a_err='integer')
    goto 9999
  end if

  diag = a%get_diag(info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_sp_getdiag')
    goto 9999
  end if

  if (do_timings) call psb_tic(idx_soc1_p0)
  call a%cp_to(acsr)
  if (do_timings) call psb_toc(idx_soc1_p0)
  if (clean_zeros) call acsr%clean_zeros(info)
  if (iorder == amg_aggr_ord_nat_) then 
    !$omp parallel do private(i)
    do i=1, nr
      ilaggr(i) = -(nr+1)
      idxs(i)   = i
      ioffs(i)  = 0
    end do
    !$omp end parallel do 
  else
    !$omp parallel do private(i)
    do i=1, nr
      ilaggr(i) = -(nr+1)
      ideg(i)   = acsr%irp(i+1) - acsr%irp(i)
      ioffs(i)  = 0
    end do
    !$omp end parallel do 
    call psb_msort(ideg,ix=idxs,dir=psb_sort_down_)
  end if
  if (do_timings) call psb_tic(idx_soc1_p1)

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

    
    !$omp parallel shared(bnds,ioffs,locnaggr,ilaggr,nr,naggr,diag,theta,nths) &
    !$omp private(icol,val,myth,kk)
    block
      integer(psb_ipk_) :: ii,nlp,k,kp,n,ia,isz, nc, i,j,m, nz, ilg,  ip, rsz
      nths = omp_get_num_threads()
      myth = omp_get_thread_num()
      rsz = nr/nths
      if (myth < mod(nr,nths)) rsz = rsz + 1
      !$omp master
      allocate(bnds(0:nths),locnaggr(0:nths+1))
      locnaggr(:) = 0
      bnds(0) = 1
      !$omp end master      
      !$omp barrier
      bnds(myth+1) = rsz
      !$omp barrier
      !$omp master
      do i=1,nths
        bnds(i) = bnds(i) + bnds(i-1)
      end do
      !$omp end master
      !$omp barrier

      !$omp  do schedule(static)
      do kk=0, nths-1
        step1: do ii=bnds(kk), bnds(kk+1)-1
          if (info /= 0) cycle
          i = idxs(ii)
          if ((i<1).or.(i>nr)) then
            info=psb_err_internal_error_
            call psb_errpush(info,name)
            cycle step1
            !goto 9999
          end if

          if (ilaggr(i) == -(nr+1)) then
            nz         = (acsr%irp(i+1)-acsr%irp(i))
            if ((nz<0).or.(nz>size(icol))) then
              info=psb_err_internal_error_
              call psb_errpush(info,name)
              cycle step1
              !goto 9999
            end if

            icol(1:nz) = acsr%ja(acsr%irp(i):acsr%irp(i+1)-1)
            val(1:nz)  = acsr%val(acsr%irp(i):acsr%irp(i+1)-1) 

            !
            ! Build the set of all strongly coupled nodes 
            !
            ip = 0
            do k=1, nz
              j   = icol(k)
              ! If any of the neighbours is already assigned,
              ! we will not reset. 
              if (ilaggr(j) > 0) cycle step1
              if (abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j)))) then
                ip = ip + 1
                icol(ip) = icol(k)
              end if
            enddo

            !
            ! If the whole strongly coupled neighborhood of I is
            ! as yet unconnected, turn it into the next aggregate.
            ! Same if ip==0 (in which case, neighborhood only
            ! contains I even if it does not look like it from matrix)
            !
            disjoint = all(ilaggr(icol(1:ip)) == -(nr+1)).or.(ip==0)
            if (disjoint) then       
              !$omp critical(update_ilaggr)
              disjoint = all(ilaggr(icol(1:ip)) == -(nr+1)).or.(ip==0)
              if (disjoint) then       
                locnaggr(kk)     = locnaggr(kk) + 1
                do k=1, ip
                  ilaggr(icol(k)) = bnds(kk)-1+locnaggr(kk)
                  ioffs(icol(k))  = kk
                end do
                ilaggr(i) = bnds(kk)-1+locnaggr(kk)
                ioffs(i)  = kk
              end if
              !$omp end critical(update_ilaggr)
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
      !$omp end master 
      !$omp barrier
      !$omp  do schedule(static) 
      do kk=0, nths-1
        do ii=bnds(kk), bnds(kk+1)-1
          if (ilaggr(ii) > 0) then
            kp = ioffs(ii) 
            ilaggr(ii) = ilaggr(ii)- (bnds(kp)-1) + locnaggr(kp)
          end if
        end do
      end do
      !$omp end do
    end block
    !$omp end parallel
  end block
#else
  step1: do ii=1, nr
    if (info /= 0) cycle
    i = idxs(ii)
    if ((i<1).or.(i>nr)) then
      info=psb_err_internal_error_
      call psb_errpush(info,name)
      cycle step1
      !goto 9999
    end if

    if (ilaggr(i) == -(nr+1)) then
      nz         = (acsr%irp(i+1)-acsr%irp(i))
      if ((nz<0).or.(nz>size(icol))) then
        info=psb_err_internal_error_
        call psb_errpush(info,name)
        cycle step1
        !goto 9999
      end if

      icol(1:nz) = acsr%ja(acsr%irp(i):acsr%irp(i+1)-1)
      val(1:nz)  = acsr%val(acsr%irp(i):acsr%irp(i+1)-1) 

      !
      ! Build the set of all strongly coupled nodes 
      !
      ip = 0
      do k=1, nz
        j   = icol(k)
        if ((1<=j).and.(j<=nr)) then 
          if (abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j)))) then
            ip = ip + 1
            icol(ip) = icol(k)
          end if
        end if
      enddo

      !
      ! If the whole strongly coupled neighborhood of I is
      ! as yet unconnected, turn it into the next aggregate.
      ! Same if ip==0 (in which case, neighborhood only
      ! contains I even if it does not look like it from matrix)
      !
      disjoint = all(ilaggr(icol(1:ip)) == -(nr+1)).or.(ip==0)
      if (disjoint) then       
        naggr     = naggr + 1
        do k=1, ip
          ilaggr(icol(k)) = naggr
        end do
        ilaggr(i) = naggr
      end if
    endif
  enddo step1
#endif
  if (debug_level >= psb_debug_outer_) then 
    write(debug_unit,*) me,' ',trim(name),&
         & ' Check   1:',naggr,count(ilaggr(1:nr) == -(nr+1)), count(ilaggr(1:nr)>0),&
         & count(ilaggr(1:nr) == -(nr+1))+count(ilaggr(1:nr)>0),nr
  end if
  if (do_timings) call psb_toc(idx_soc1_p1)
  if (do_timings) call psb_tic(idx_soc1_p2)
  !
  ! Phase two: join the neighbours
  !
  !$omp workshare
  tmpaggr = ilaggr
  !$omp end workshare
  !$omp parallel do schedule(static) shared(tmpaggr,ilaggr,nr,naggr,diag,theta)& 
  !$omp     private(ii,i,j,k,nz,icol,val,ip)
  step2: do ii=1,nr
    i = idxs(ii)

    if (ilaggr(i) == -(nr+1)) then         
      nz         = (acsr%irp(i+1)-acsr%irp(i))
      if (nz == 1) cycle step2
      icol(1:nz) = acsr%ja(acsr%irp(i):acsr%irp(i+1)-1)
      val(1:nz)  = acsr%val(acsr%irp(i):acsr%irp(i+1)-1) 

      !
      ! Find the most strongly connected neighbour that is
      ! already aggregated, if any, and join its aggregate
      !
      cpling = dzero
      ip = 0
      do k=1, nz
        j   = icol(k)
        if ((1<=j).and.(j<=nr)) then 
          if ((abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j))))&
               & .and. (tmpaggr(j) > 0).and. (abs(val(k)) > cpling)) then
            ip = k
            cpling = abs(val(k))
          end if
        end if
      enddo
      if (ip > 0) then 
        ilaggr(i) = ilaggr(icol(ip))
      end if
    end if
  end do step2
  !$omp end parallel do
  if (do_timings) call psb_toc(idx_soc1_p2)
  if (debug_level >= psb_debug_outer_) then 
    write(debug_unit,*) me,' ',trim(name),&
         & ' Check 1.5:',naggr,count(ilaggr(1:nr) == -(nr+1)), count(ilaggr(1:nr)>0),&
         & count(ilaggr(1:nr) == -(nr+1))+count(ilaggr(1:nr)>0),nr
  end if

  if (do_timings) call psb_tic(idx_soc1_p3)
  !
  ! Phase three: sweep over leftovers, if any 
  !
  step3: do ii=1,nr
    i = idxs(ii)

    if (ilaggr(i) < 0) then
      nz         = (acsr%irp(i+1)-acsr%irp(i))
      if (nz == 1) cycle step3
      icol(1:nz) = acsr%ja(acsr%irp(i):acsr%irp(i+1)-1)
      val(1:nz)  = acsr%val(acsr%irp(i):acsr%irp(i+1)-1) 
      !
      ! Find its strongly  connected neighbourhood not 
      ! already aggregated, and make it into a new aggregate.
      !
      cpling = dzero
      ip = 0
      do k=1, nz
        j   = icol(k)
        if ((1<=j).and.(j<=nr)) then 
          if ((abs(val(k)) > theta*sqrt(abs(diag(i)*diag(j))))&
               & .and. (ilaggr(j) < 0))  then
            ip = ip + 1
            icol(ip) = icol(k)
          end if
        end if
      enddo
      if (ip > 0) then
        naggr     = naggr + 1
        ilaggr(i) = naggr
        do k=1, ip
          ilaggr(icol(k)) = naggr
        end do
      else
        !
        ! This should not happen: we did not even connect with ourselves,
        ! but it's not a singleton. 
        !
        naggr     = naggr + 1
        ilaggr(i) = naggr
      end if
    end if
  end do step3

  ! Any leftovers?
  do i=1, nr
    if (ilaggr(i) < 0) then
      nz = (acsr%irp(i+1)-acsr%irp(i))
      if (nz == 1) then
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
  if (do_timings) call psb_toc(idx_soc1_p3)
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
  if (debug_level >= psb_debug_outer_) then 
    write(debug_unit,*) me,' ',trim(name),&
         & ' Check   2:',naggr,count(ilaggr(1:nr) == -(nr+1)), count(ilaggr(1:nr)>0),&
         & count(ilaggr(1:nr) == -(nr+1))+count(ilaggr(1:nr)>0),nr
  end if

  call acsr%free()
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_d_soc1_map_bld


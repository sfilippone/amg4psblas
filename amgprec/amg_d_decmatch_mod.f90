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
module amg_d_decmatch_mod

  use iso_c_binding
  use psb_base_cbind_mod

  interface new_Match_If
    function dnew_Match_If(ipar,matching,lambda,nr, irp, ja, val, diag, w, mate) &
         & bind(c,name="dnew_Match_If") result(res)
      use iso_c_binding
      import :: psb_c_ipk_, psb_c_lpk_, psb_c_mpk_, psb_c_epk_
      implicit none

      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nr,ipar,matching
      real(c_double), value :: lambda      
      type(c_ptr), value :: irp, ja, mate
      type(c_ptr), value :: val, diag, w 
    end function dnew_Match_If
  end interface new_Match_If

  interface amg_build_decmatch
    module procedure amg_dbuild_decmatch
  end interface amg_build_decmatch

  logical, parameter, private :: print_statistics=.false.
contains

  subroutine amg_ddecmatch_build_prol(w,a,desc_a,ilaggr,nlaggr,prol,info,&
       & symmetrize,reproducible,display_inp, display_out, print_out, &
       & parallel, matching,lambda)
    use psb_base_mod
    use psb_util_mod
    use iso_c_binding
    implicit none
    real(psb_dpk_), allocatable, intent(inout)  :: w(:)
    type(psb_dspmat_type), intent(inout)        :: a
    type(psb_desc_type)                         :: desc_a
    integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:)
    integer(psb_lpk_), allocatable, intent(out) :: nlaggr(:)
    type(psb_ldspmat_type), intent(out)         :: prol
    integer(psb_ipk_), intent(out)              :: info
    logical, optional, intent(in)               :: display_inp, display_out, reproducible
    logical, optional, intent(in)               :: symmetrize, print_out, parallel
    integer(psb_ipk_), optional, intent(in)     :: matching
    real(psb_dpk_), optional, intent(in)        :: lambda

    !
    !
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: iam, np, iown
    integer(psb_ipk_) :: nr, nc, sweep, nzl, ncsave, nct, idx
    integer(psb_lpk_) :: i, k, kg, idxg, ntaggr, naggrm1, naggrp1, &
         & ip, nlpairs, nlsingl, nunmatched, lnr
    real(psb_dpk_)    :: wk, widx, wmax, nrmagg
    real(psb_dpk_), allocatable    :: wtemp(:)
    integer(psb_ipk_), allocatable :: mate(:)
    integer(psb_lpk_), allocatable ::  ilv(:)
    integer(psb_ipk_), save     :: cnt=1
    character(len=256)          :: aname
    type(psb_ld_coo_sparse_mat) :: tmpcoo
    logical  :: display_out_, print_out_, reproducible_, parallel_
    integer(psb_ipk_) :: matching_
    real(psb_dpk_)    :: lambda_    
    logical, parameter  :: dump=.false., debug=.false., dump_mate=.false., &
         & debug_ilaggr=.false., debug_sync=.false.
    integer(psb_ipk_), save :: idx_bldmtc=-1, idx_phase1=-1, idx_phase2=-1, idx_phase3=-1
    logical, parameter :: do_timings=.true.

    ictxt = desc_a%get_ctxt()
    call psb_info(ictxt,iam,np)

    if ((do_timings).and.(idx_phase1==-1))       &
         & idx_phase1 = psb_get_timer_idx("MBP_BLDP: phase1       ")
    if ((do_timings).and.(idx_bldmtc==-1))       &
         & idx_bldmtc = psb_get_timer_idx("MBP_BLDP: buil_matching")
    if ((do_timings).and.(idx_phase2==-1))       &
         & idx_phase2 = psb_get_timer_idx("MBP_BLDP: phase2       ")
    if ((do_timings).and.(idx_phase3==-1))       &
         & idx_phase3 = psb_get_timer_idx("MBP_BLDP: phase3       ")

    if (do_timings) call psb_tic(idx_phase1)

    if (present(display_out)) then
      display_out_ = display_out
    else
      display_out_ = .false.
    end if
    if (present(print_out)) then
      print_out_ = print_out
    else
      print_out_ = .false.
    end if
    if (present(reproducible)) then
      reproducible_ = reproducible
    else
      reproducible_ = .false.
    end if

    if (present(parallel)) then
      parallel_ = parallel
    else
      parallel_ = .true.
    end if

    if (present(matching)) then
      matching_ = matching
    else
      matching_ = 2
    end if

    if (present(lambda)) then
      lambda_ = lambda
    else
      lambda_ = 2.0
    end if

    allocate(nlaggr(0:np-1),stat=info)
    if (info /= 0) then
      return
    end if

    nlaggr = 0
    ilv    = [(i,i=1,desc_a%get_local_cols())]
    call desc_a%l2gip(ilv,info,owned=.false.)

    call psb_geall(ilaggr,desc_a,info)
    ilaggr = -1
    call psb_geasb(ilaggr,desc_a,info)
    nr = a%get_nrows()
    nc = a%get_ncols()
    if (size(w) < nc) then
      call psb_realloc(nc,w,info)
    end if
    call psb_halo(w,desc_a,info)

    if (debug) write(0,*) iam,' buildprol into buildmatching:',&
         & nr, nc
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (iam == 0) write(0,*)' buildprol into buildmatching:',&
           & nr, nc
    end if
    if (do_timings) call psb_toc(idx_phase1)
    if (do_timings) call psb_tic(idx_bldmtc)
    call amg_dbuild_decmatch(parallel_,matching_,lambda_,w,a,desc_a,mate,info)
    if (do_timings) call psb_toc(idx_bldmtc)
    if (debug) write(0,*) iam,' buildprol from buildmatching:',&
         & info
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (iam == 0) write(0,*)' out from buildmatching:', info
    end if

    if (info == 0) then
      if (do_timings) call psb_tic(idx_phase2)
      if (debug_sync) then
        call psb_barrier(ictxt)
        if (iam == 0) write(0,*)' Into building the tentative prol:'
      end if

      call psb_geall(wtemp,desc_a,info)
      wtemp      = dzero
      call psb_geasb(wtemp,desc_a,info)

      nlaggr(iam) = 0
      nlpairs    = 0
      nlsingl    = 0
      nunmatched = 0
      !
      ! First sweep
      ! On return from build_matching, mate has been converted to local numbering,
      ! so assigning to idx is OK.
      !
      do k=1, nr
        idx = mate(k)
        !
        ! Figure out an allocation of aggregates to processes
        !
        if (idx < 0) then
          !
          !  Unmatched vertex, potential  singleton.
          !
          nunmatched = nunmatched + 1
          if (abs(w(k))<epsilon(nrmagg)) then
            ! Keep it unaggregated
            wtemp(k) = dzero
          else
            ! Create a singleton aggregate
            nlaggr(iam) = nlaggr(iam) + 1
            ilaggr(k)   = nlaggr(iam)
            wtemp(k)    = w(k)/abs(w(k))
            nlsingl     = nlsingl + 1
          end if
!!$          write(0,*) k,mate(k),ilaggr(k),' negative match ',abs(w(k)), epsilon(nrmagg)
        else  if (idx > nc) then
          write(0,*) 'Impossible: mate(k) > nc'
          cycle
        else

          if (ilaggr(k) == -1) then

            wk   = w(k)
            widx = w(idx)
            wmax = max(abs(wk),abs(widx))
            nrmagg = wmax*sqrt((wk/wmax)**2+(widx/wmax)**2)
            if (nrmagg > epsilon(nrmagg)) then
              if (idx <= nr) then
                if (ilaggr(idx) == -1) then
                  ! Now, if both vertices are local, the aggregate is local
                  ! (kinda obvious).
                  nlaggr(iam) = nlaggr(iam) + 1
                  ilaggr(k)   = nlaggr(iam)
                  ilaggr(idx) = nlaggr(iam)
                  wtemp(k)    = w(k)/nrmagg
                  wtemp(idx)  = w(idx)/nrmagg
                end if
                nlpairs = nlpairs+1
              else
                write(0,*) 'Really?  mate(k) > nr? ',mate(k),nr
              end if
            else
              if (abs(w(k))<epsilon(nrmagg)) then
                ! Keep it unaggregated
                wtemp(k) = dzero
              else
                ! Create a singleton aggregate
                nlaggr(iam) = nlaggr(iam) + 1
                ilaggr(k)   = nlaggr(iam)
                wtemp(k)    = w(k)/abs(w(k))
                nlsingl     = nlsingl + 1
              end if
            end if
          end if
        end if
      end do
      if (do_timings) call psb_toc(idx_phase2)
      if (do_timings) call psb_tic(idx_phase3)

      ! Ok, now compute offsets, gather halo and fix non-local
      ! aggregates  (those where ilaggr == -2)
      call psb_sum(ictxt,nlaggr)
      ntaggr  = sum(nlaggr(0:np-1))
      naggrm1 = sum(nlaggr(0:iam-1))
      naggrp1 = sum(nlaggr(0:iam))
      !
      ! Shift all indices already assigned (i.e. >0)
      !
      do k=1,nr
        if (ilaggr(k) > 0) then
          ilaggr(k) = ilaggr(k) + naggrm1
!!$        else
!!$          write(0,*) 'Leftover ILAGGR',k,ilaggr(k),mate(k),abs(w(k)),epsilon(nrmagg)
        end if
      end do
      call psb_halo(ilaggr,desc_a,info)
      call psb_halo(wtemp,desc_a,info)
      ! Cleanup as yet unmarked entries
      do k=1,nr
        if (ilaggr(k) == -2) then
          idx = mate(k)
          if (idx > nr) then
            i   = ilaggr(idx)
            if (i > 0) then
              ilaggr(k) = i
            else
              write(0,*) 'Error : unresolved (paired) index ',k,idx,i,nr,nc, ilv(k),ilv(idx)
            end if
          else
            write(0,*) 'Error : unresolved (paired) index ',k,idx,i,nr,nc, ilv(k),ilv(idx)
          end if
        end if
        if (ilaggr(k) <0) then
          write(0,*) 'Decmatch: Funny number: ',k,ilv(k),ilaggr(k),wtemp(k)
        end if
      end do
      if (debug_sync) then
        call psb_barrier(ictxt)
        if (iam == 0) write(0,*)' Done building the tentative prol:'
      end if


      if (dump_mate) then
        block
          integer(psb_lpk_), allocatable :: glaggr(:)
          write(aname,'(a,i3.3,a,i3.3,a)') 'mateg-',cnt,'-p',iam,'.mtx'
          open(20,file=aname)
          write(20,'(a,I3,a)') '% sparse vector on process ',iam,' '
          do k=1, nr
            write(20,'(3(I8,1X))') ilv(k),ilv(mate(k))
          end do
          close(20)
          write(aname,'(a,i3.3,a,i3.3,a)') 'nloc-',cnt,'-p',iam,'.mtx'
          open(20,file=aname)
          write(20,'(a,I3,a)') '% sparse vector on process ',iam,' '
          write(20,'(a,I12,a)') 'nlpairs ',nlpairs
          write(20,'(a,I12,a)') 'nlsingl ',nlsingl
          write(20,'(a,I12,a)') 'nlaggr(iam) ',nlaggr(iam)
          close(20)
          write(aname,'(a,i3.3,a,i3.3,a,i3.3,a)') 'ilaggr-',cnt,'-i',iam,'-p',np,'.mtx'
          open(20,file=aname)
          write(20,'(a,I3,a)') '% sparse vector on process ',iam,' '
          do k=1, nr
            write(20,'(3(I8,1X))') ilv(k),ilaggr(k)
          end do
          close(20)

          write(aname,'(a,i3.3,a,i3.3,a)') 'glaggr-',cnt,'-p',np,'.mtx'
          call psb_gather(glaggr,ilaggr,desc_a,info,root=izero)
          if (iam==0) call mm_array_write(glaggr,'Aggregates ',info,filename=aname)

          cnt=cnt+1
        end block
      end if
      block
        integer(psb_lpk_) :: v(3)
        v(1) = nunmatched
        v(2) = nlsingl
        v(3) = nlpairs
        call psb_sum(ictxt,v)
        nunmatched = v(1)
        nlsingl    = v(2)
        nlpairs    = v(3)

      end block
      if (print_statistics) then
        if (iam == 0) then
          write(0,*) 'Matching statistics: Unmatched nodes ',&
               & nunmatched,' Singletons:',nlsingl,' Pairs:',nlpairs
        end if
      end if
      
      if (display_out_) then
        block
          integer(psb_ipk_) :: idx
          !
          !  And finally print out
          !
          do i=0,np-1
            call psb_barrier(ictxt)
            if (iam == i) then
              write(0,*) 'Process ', iam,' hosts aggregates:  (',naggrm1+1,' : ',naggrp1-1,')'
              do k=1, nr
                idx  = mate(k)
                kg    = ilv(k)
                if (idx >0) then
                  idxg  = ilv(idx)
                else
                  idxg = -1
                end if
                if (idx < 0) then
                  write(0,*) kg,': singleton   (',kg,' ( Proc',iam,') ) into aggregate => ', ilaggr(k)
                else if (idx <= nr) then
                  write(0,*) kg,': paired with (',idxg,' ( Proc',iam,') ) into aggregate => ', ilaggr(k)
                else
                  call desc_a%indxmap%qry_halo_owner(idx,iown,info)
                  write(0,*) kg,': paired with (',idxg,' ( Proc',iown,') ) into aggregate => ', ilaggr(k)
                end if
              end do
              flush(0)
            end if
          end do
        end block
      end if

      ! Dirty trick: allocate tmpcoo with local
      ! number of aggregates, then change to ntaggr.
      ! Just to make sure the allocation is not global
      lnr = nr
      call tmpcoo%allocate(lnr,nlaggr(iam),lnr)
      k = 0
      do i=1,nr
        !
        ! Note: at this point, a value ilaggr(i)<=0
        ! tags an unaggregated row, and it has to be
        ! left alone (i.e.: it should stay at fine level only)
        !
        if (ilaggr(i)>0) then
          k = k + 1
          tmpcoo%val(k) = wtemp(i)
          tmpcoo%ia(k)  = i
          tmpcoo%ja(k)  = ilaggr(i)
        end if
      end do
      call tmpcoo%set_nzeros(k)
      call tmpcoo%set_dupl(psb_dupl_add_)
      call tmpcoo%set_sorted() ! This is now in row-major

      if (display_out_) then
        call psb_barrier(ictxt)
        flush(0)

        if (iam == 0) write(0,*) 'Prolongator: '
        flush(0)
        do i=0,np-1
          call psb_barrier(ictxt)
          if (iam == i) then
            do k=1, nr
              write(0,*) ilv(tmpcoo%ia(k)),tmpcoo%ja(k), tmpcoo%val(k)
            end do
            flush(0)
          end if
        end do
      end if

      call prol%mv_from(tmpcoo)
      if (do_timings) call psb_toc(idx_phase3)

      if (print_out_) then
        write(aname,'(a,i3.3,a)') 'prol-g-',iam,'.mtx'
        call prol%print(fname=aname,head='Test ',ivr=ilv)
        write(aname,'(a,i3.3,a)') 'prol-',iam,'.mtx'
        call prol%print(fname=aname,head='Test ')
      end if

    else
      write(0,*) iam,' :  error from Matching: ',info
    end if

  end subroutine amg_ddecmatch_build_prol

  subroutine amg_dbuild_decmatch(parallel,matching,lambda,w,a,desc_a,mate,info)
    use psb_base_mod
    use psb_util_mod
    use iso_c_binding
    implicit none
    logical, intent(in)               :: parallel
    integer(psb_ipk_), intent(in)     :: matching
    real(psb_dpk_), intent(in)        :: lambda
    real(psb_dpk_), target                      :: w(:)
    type(psb_dspmat_type), intent(inout)        :: a
    type(psb_desc_type)                         :: desc_a
    integer(psb_ipk_), allocatable, intent(out), target :: mate(:)
    integer(psb_ipk_), intent(out)              :: info

    type(psb_d_csr_sparse_mat), target :: tcsr
    real(psb_dpk_), allocatable, target :: diag(:)
    real(psb_dpk_) :: ph0t,ph1t,ph2t
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: iam, np
    integer(psb_ipk_)  :: nr, nc, nz, i, nunmatch, ipar 
    integer(psb_ipk_), save :: cnt=2
    logical, parameter :: debug=.false., dump_ahat=.false., debug_sync=.false.
    logical, parameter :: old_style=.false., sort_minp=.true.
    character(len=40) :: name='build_matching', fname
    integer(psb_ipk_), save :: idx_cmboxp=-1, idx_bldahat=-1, idx_phase2=-1, idx_phase3=-1
    logical, parameter :: do_timings=.true.

    ictxt = desc_a%get_ctxt()
    call psb_info(ictxt,iam,np)

    nr = a%get_nrows()
    call a%cp_to(tcsr)
    call psb_realloc(nr,mate,info)
    diag = a%get_diag(info)
    if (parallel) then
      ipar = 2
    else
      ipar = 1
    end if
    !
    ! Now call matching!
    !
    if (debug) write(0,*) iam,' buildmatching into PMatchBox:'
    if (do_timings) call psb_tic(idx_cmboxp)
    info = dnew_Match_If(ipar,matching,lambda,nr,c_loc(tcsr%irp),c_loc(tcsr%ja),&
         & c_loc(tcsr%val),c_loc(diag),c_loc(w),c_loc(mate))
    if (do_timings) call psb_toc(idx_cmboxp)
    if (debug) write(0,*) iam,' buildmatching from PMatchBox:', info
    if (debug_sync) then
      call psb_max(ictxt,info)
      if (iam == 0) write(0,*)' done PMatchBox', info
    end if

    if (do_timings) call psb_tic(idx_phase3)
    nunmatch = count(mate(1:nr)<=0)
        ! call psb_sum(ictxt,nunmatch)
    if (nunmatch /= 0) write(0,*) iam,' Unmatched nodes local imbalance ',nunmatch
 !    if (count(mate(1:nr)<0) /= nunmatch) write(0,*) iam,' Matching results ?',&
 !        & nunmatch, count(mate(1:nr)<0)
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (iam == 0) write(0,*)' done build_matching '
    end if

    if (do_timings) call psb_toc(idx_phase3)
    return

9999 continue
    call psb_error(ictxt)

  end subroutine amg_dbuild_decmatch

end module amg_d_decmatch_mod

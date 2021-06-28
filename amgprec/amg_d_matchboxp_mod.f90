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
! moved here from amg4psblas-extension
!
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
module dmatchboxp_mod

  use iso_c_binding
  use psb_base_cbind_mod

  interface MatchBoxPC
    subroutine dMatchBoxPC(nlver,nledge,verlocptr,verlocind,edgelocweight,&
         & verdistance, mate, myrank, numprocs, icomm,&
         & msgindsent,msgactualsent,msgpercent,&
         & ph0_time, ph1_time, ph2_time, ph1_card, ph2_card) bind(c,name='dMatchBoxPC')
      use iso_c_binding
      import :: psb_c_ipk_, psb_c_lpk_, psb_c_mpk_, psb_c_epk_
      implicit none

      integer(psb_c_lpk_), value :: nlver,nledge
      integer(psb_c_mpk_), value :: myrank, numprocs, icomm
      integer(psb_c_lpk_) :: verlocptr(*),verlocind(*), verdistance(*)
      integer(psb_c_lpk_) :: mate(*)
      integer(psb_c_lpk_) :: msgindsent(*),msgactualsent(*)
      real(c_double)      :: ph0_time, ph1_time, ph2_time
      integer(psb_c_lpk_) :: ph1_card(*),ph2_card(*)
      real(c_double) :: edgelocweight(*)
      real(c_double)     :: msgpercent(*)
    end subroutine dMatchBoxPC
  end interface MatchBoxPC

  interface i_aggr_assign
    module procedure i_daggr_assign
  end interface i_aggr_assign

  interface build_matching
    module procedure dbuild_matching
  end interface build_matching

  interface build_ahat
    module procedure dbuild_ahat
  end interface build_ahat

  interface psb_gtranspose
    module procedure psb_dgtranspose
  end interface psb_gtranspose

  interface psb_htranspose
    module procedure psb_dhtranspose
  end interface psb_htranspose

  interface PMatchBox
    module procedure dPMatchBox
  end interface PMatchBox

  logical, parameter, private :: print_statistics=.false.
contains

  subroutine dmatchboxp_build_prol(w,a,desc_a,ilaggr,nlaggr,prol,info,&
       & symmetrize,reproducible,display_inp, display_out, print_out)
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
    logical, optional, intent(in)               :: symmetrize, print_out

    !
    !
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: iam, np, iown
    integer(psb_ipk_) :: nr, nc, sweep, nzl, ncsave, nct, idx
    integer(psb_lpk_) :: i, k, kg, idxg, ntaggr, naggrm1, naggrp1, &
         & ip, nlpairs, nlsingl, nunmatched, lnr
    real(psb_dpk_)    :: wk, widx, wmax, nrmagg
    real(psb_dpk_), allocatable    :: wtemp(:)
    integer(psb_lpk_), allocatable :: mate(:), ilv(:)
    integer(psb_ipk_), save     :: cnt=1
    character(len=256)          :: aname
    type(psb_ld_coo_sparse_mat) :: tmpcoo
    logical  :: display_out_, print_out_, reproducible_
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
    call build_matching(w,a,desc_a,mate,info,display_inp=display_inp,symmetrize=symmetrize)
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
              else if (idx <= nc) then
                !
                ! This pair involves a non-local vertex.
                ! Set wtemp, then apply a tie-breaking algorithm
                wtemp(k)    = w(k)/nrmagg
                idxg  = ilv(idx)
                kg    = ilv(k)
                if (reproducible_) then
                  !
                  ! Tie-break by assigning to the lowest-index process
                  ! so that the numbering is the same as with NP=1
                  !
                  if (kg < idxg) then
                    nlaggr(iam) = nlaggr(iam) + 1
                    ilaggr(k)   = nlaggr(iam)
                    nlpairs = nlpairs+1
                  else
                    ilaggr(k) = -2
                  end if
                else
                  ! Use a statistically unbiased tie-breaking rule,
                  ! this will give an even spread.
                  ! Delegate to i_aggr_assign.
                  ! Should be a symmetric function.
                  !
                  call desc_a%indxmap%qry_halo_owner(idx,iown,info)
                  ip    = i_aggr_assign(iam, iown, kg, idxg, wk, widx, nrmagg)
                  if (iam == ip) then
                    nlaggr(iam) = nlaggr(iam) + 1
                    ilaggr(k)   = nlaggr(iam)
                    nlpairs = nlpairs+1
                  else
                    ilaggr(k) = -2
                  end if
                end if
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
        if (ilaggr(k) > 0) ilaggr(k) = ilaggr(k) + naggrm1
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
          write(0,*) 'Matchboxp: Funny number: ',k,ilv(k),ilaggr(k),wtemp(k)
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

  end subroutine dmatchboxp_build_prol

  function i_daggr_assign(iam, iown, kg, idxg, wk, widx, nrmagg) &
       & result(iproc)
    !
    ! How to break ties? This
    ! must  be a symmetric function, i.e.
    ! the result value iproc MUST be the same upon
    ! swapping simultaneously all  pairs
    ! (iam,iown) (kg,idxg) (wk,widx)
    !
    implicit none
    integer(psb_ipk_)             :: iproc
    integer(psb_ipk_), intent(in) :: iam, iown
    integer(psb_lpk_), intent(in) :: kg, idxg
    real(psb_dpk_), intent(in)    :: wk, widx, nrmagg
    !
    integer(psb_lpk_) :: kg2, idxg2

    idxg2 = mod(idxg,2)
    kg2   = mod(kg,2)

    !
    ! In this particular tie-breaking rule we are ignoring WK, WIDX.
    ! This should statistically entail an even spread of aggregates.
    !
    if ((kg2/=0).and.(idxg2/=0)) then
      ! If both  global indices are  odd,
      !    assign to the higher number process
      iproc = max(iam,iown)

    else if ((kg2==0).and.(idxg2==0)) then
      ! If both global indices are even,
      !    assign to the lower number process;
      iproc = min(iam,iown)
    else
      ! If the global indices  are mixed,
      !    then  assign to the owner of the odd index
      if (kg2 /= 0) then
        iproc = iam
      else
        iproc = iown
      end if
    end if
  end function i_daggr_assign


  subroutine dbuild_matching(w,a,desc_a,mate,info,display_inp, symmetrize)
    use psb_base_mod
    use psb_util_mod
    use iso_c_binding
    implicit none
    real(psb_dpk_)                              :: w(:)
    type(psb_dspmat_type), intent(inout)       :: a
    type(psb_desc_type)                         :: desc_a
    integer(psb_lpk_), allocatable, intent(out) :: mate(:)
    integer(psb_ipk_), intent(out)              :: info
    logical, optional                           :: display_inp, symmetrize

    type(psb_dspmat_type) :: ahatnd
    type(psb_ld_csr_sparse_mat) :: tcsr
    type(psb_desc_type) :: desc_blk
    integer(psb_lpk_), allocatable :: vnl(:)
    integer(psb_lpk_), allocatable :: ph1crd(:), ph2crd(:)
    integer(psb_lpk_), allocatable :: msgis(:),msgas(:)
    real(psb_dpk_), allocatable :: msgprc(:)
    real(psb_dpk_) :: ph0t,ph1t,ph2t
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: iam, np
    integer(psb_lpk_) :: nr, nc, nz, i, nunmatch
    integer(psb_ipk_), save :: cnt=2
    logical, parameter :: debug=.false., dump_ahat=.false., debug_sync=.false.
    logical, parameter :: old_style=.false., sort_minp=.true.
    character(len=40) :: name='build_matching', fname
    integer(psb_ipk_), save :: idx_cmboxp=-1, idx_bldahat=-1, idx_phase2=-1, idx_phase3=-1
    logical, parameter :: do_timings=.true.

    ictxt = desc_a%get_ctxt()
    call psb_info(ictxt,iam,np)

    if (debug) write(0,*) iam,' buildmatching into build_ahat:'
    if ((do_timings).and.(idx_bldahat==-1))       &
         & idx_bldahat = psb_get_timer_idx("BLD_MTCH: bld_ahat")
    if ((do_timings).and.(idx_cmboxp==-1))       &
         & idx_cmboxp = psb_get_timer_idx("BLD_MTCH: PMatchBox")
    if ((do_timings).and.(idx_phase2==-1))       &
         & idx_phase2 = psb_get_timer_idx("BLD_MTCH: phase2")
    if ((do_timings).and.(idx_phase3==-1))       &
         & idx_phase3 = psb_get_timer_idx("BLD_MTCH: phase3")

    if (debug) write(0,*) iam,' buildmatching from cd_renum:',info
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (iam == 0) write(0,*)' Into build_ahat:'
    end if
    if (do_timings) call psb_tic(idx_bldahat)
    call build_ahat(w,a,ahatnd,desc_a,info,symmetrize=symmetrize)
    if (do_timings) call psb_toc(idx_bldahat)
    if (info /= 0) then
      write(0,*) 'Error from build_ahat ', info
    end if
    if (debug) write(0,*) iam,' buildmatching from build_ahat:',info
    if (dump_ahat) then
      block
        character(len=40) :: fname
        integer(psb_ipk_) :: k, nr
        type(psb_ldspmat_type) :: tmp_mat
        integer(psb_lpk_), allocatable :: ilv(:)
        if (.false.) then
          ilv = desc_a%get_global_indices(owned=.false.)
          write(fname,'(a,i3.3,a,i3.3,a)') 'w-i',cnt,'-p',iam,'.mtx'
          open(20,file=fname)
          write(20,'(a,I3,a)') '% W vector           ',iam,' '
          do k=1, nr
            write(20,'(I8,1X,es26.18)') ilv(k),w(k)
          end do
          close(20)
          write(fname,'(a,i3.3,a,i3.3,a)') 'aa-i',cnt,'-p',iam,'.mtx'
!!$          call a%print(fname=fname,head='Original matrix  ',iv=ilv)
          write(fname,'(a,i3.3,a,i3.3,a)') 'ahat-i',cnt,'-p',iam,'.mtx'
          call ahatnd%print(fname=fname,head='Input to matching ')
        else
          write(fname,'(a,i3.3,a,i3.3,a)') 'ahat-i',cnt,'-p',np,'.mtx'
          call psb_gather(tmp_mat,ahatnd,desc_a,info,root=izero)
          if (iam == 0) call tmp_mat%print(fname=fname,head='Input to matching ')
        end if

      end block
    end if
    if (do_timings) call psb_tic(idx_phase2)
    !
    ! Now AHATND should be symmetric and positive, without a diagonal.
    ! Almost ready to call matching routine.
    !
    call ahatnd%mv_to(tcsr)
    nr = tcsr%get_nrows()
    nc = tcsr%get_ncols()
    nz = tcsr%get_nzeros()
    allocate(vnl(0:np),mate(max(nr,nc)),&
         & msgis(np),msgas(np),msgprc(np),&
         & ph1crd(np),ph2crd(np),stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif
    vnl = 0
    vnl(iam) = nr
    call psb_sum(ictxt,vnl)
    vnl(1:np) = vnl(0:np-1)
    vnl(0) = 1
    do i=1,np
      vnl(i) = vnl(i-1)+vnl(i)
    end do
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (iam == 0) write(0,*)' renum_blk'
    end if

    call psb_cd_renum_block(desc_a,desc_blk,info)
    associate (vlptr => tcsr%irp, vlind => tcsr%ja, ewght => tcsr%val)
      ! Put column indices in global numbering
      call psb_loc_to_glob(vlind,desc_blk,info,iact='E')
      mate          = -1

      if (sort_minp) then
        block
          integer(psb_ipk_) :: ir1,ir2,nrz
          do i=1,nr
            ir1 = vlptr(i)
            ir2 = vlptr(i+1) -1
            nrz = (ir2-ir1+1)
            call fix_order(nrz,vlind(ir1:ir2),ewght(ir1:ir2),info)
          end do
        end block
      end if

      if (debug_sync) then
        call psb_barrier(ictxt)
        if (iam == 0) write(0,*)' into PMatchBox'
      end if
      if (do_timings) call psb_toc(idx_phase2)
      !
      ! Now call matching!
      !
      if (debug) write(0,*) iam,' buildmatching into PMatchBox:'
      if (do_timings) call psb_tic(idx_cmboxp)
      call PMatchBox(nr,nz,vlptr,vlind,ewght,&
           & vnl, mate, iam, np,ictxt,&
           & msgis,msgas,msgprc,ph0t,ph1t,ph2t,ph1crd,ph2crd,info,display_inp)
      if (do_timings) call psb_toc(idx_cmboxp)
      if (debug) write(0,*) iam,' buildmatching from PMatchBox:', info
      if (debug_sync) then
        call psb_max(ictxt,info)
        if (iam == 0) write(0,*)' done PMatchBox', info
      end if
    end associate
    if (do_timings) call psb_tic(idx_phase3)
    nunmatch = count(mate(1:nr)<=0)
    call psb_glob_to_loc(mate(1:nr),desc_blk,info,iact='I',owned=.false.)
    nunmatch = abs(nunmatch - count(mate(1:nr)<=0))
    ! call psb_sum(ictxt,nunmatch)
    if (nunmatch /= 0) write(0,*) iam,' Unmatched nodes local imbalance ',nunmatch
!!$    if (count(mate(1:nr)<0) /= nunmatch) write(0,*) iam,' Matching results ?',&
!!$         & nunmatch, count(mate(1:nr)<0)
    call desc_blk%free(info)
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (iam == 0) write(0,*)' done build_matching '
    end if

    if (dump_ahat) then
      block
        integer(psb_lpk_), allocatable :: mg(:), ml(:)
        real(psb_dpk_), allocatable    :: wl(:), wg(:)
        ml=mate
        call psb_loc_to_glob(ml(1:nr),desc_a,info,iact='E')
        write(fname,'(a,i3.3,a,i3.3,a)') 'mate-i',cnt,'-p',np,'.mtx'
        call psb_gather(mg,ml,desc_a,info,root=izero)
        if (iam==0) call mm_array_write(mg,'Output from matching ',info,filename=fname)
        wl=w
        write(fname,'(a,i3.3,a,i3.3,a)') 'w-',cnt,'-p',np,'.mtx'
        call psb_gather(wg,wl,desc_a,info,root=izero)
        if (iam==0) call mm_array_write(wg,'Input smooth vector ',info,filename=fname)

      end block
      cnt = cnt + 1
    end if
    if (do_timings) call psb_toc(idx_phase3)
    return

9999 continue
    call psb_error(ictxt)
  contains
    subroutine fix_order(n,ja,val,iret)
      use psb_base_mod
      implicit none
      integer(psb_ipk_), intent(in)    :: n
      integer(psb_lpk_), intent(inout) :: ja(:)
      real(psb_dpk_), intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)    :: iret
      integer(psb_lpk_), allocatable   :: ix(:)
      real(psb_dpk_), allocatable      :: tmp(:)

      allocate(ix(n), tmp(n),stat=iret)
      if (iret /= 0) return
      call psb_msort(ja(1:n),ix=ix,dir=psb_sort_up_)
      tmp(1:n) = val(ix(1:n))
      val(1:n) = tmp(1:n)
    end subroutine fix_order

  end subroutine dbuild_matching

  subroutine dbuild_ahat(w,a,ahat,desc_a,info,symmetrize)
    use psb_base_mod
    implicit none
    real(psb_dpk_), intent(in)           :: w(:)
    type(psb_dspmat_type), intent(inout) :: a
    type(psb_dspmat_type), intent(out)   :: ahat
    type(psb_desc_type)                   :: desc_a
    integer(psb_ipk_), intent(out)        :: info
    logical, optional                     :: symmetrize

    type(psb_dspmat_type)      :: atnd
    type(psb_d_coo_sparse_mat) :: tcoo1, tcoo2, tcoo3
    real(psb_dpk_), allocatable :: diag(:)
    integer(psb_lpk_), allocatable :: ilv(:)
    logical, parameter :: debug=.false., dump=.false., dump_ahat=.false., dump_inp=.false.
    logical :: symmetrize_
    real(psb_dpk_)    :: aii, ajj, aij, wii, wjj, tmp1, tmp2, minabs, edgnrm
    integer(psb_ipk_) :: nr, nc, nz, i, nz2, nrg, ii, jj, k, k2, ncols
    integer(psb_lpk_) :: nzglob
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: me, np
    character(len=80) :: aname
    real(psb_dpk_), parameter :: eps=epsilon(1.d0)
    integer(psb_ipk_), save   :: idx_glbt=-1, idx_phase1=-1, idx_phase2=-1
    logical, parameter :: do_timings=.true.
    logical, parameter :: debug_symmetry = .false., check_size=.false.
    logical, parameter :: unroll_logtrans=.false.

    ictxt = desc_a%get_ctxt()
    call psb_info(ictxt,me,np)
    if (present(symmetrize)) then
      symmetrize_ = symmetrize
    else
      symmetrize_ = .false.
    end if
    if ((do_timings).and.(idx_phase1==-1))       &
         & idx_phase1 = psb_get_timer_idx("BLD_AHAT: phase1       ")
    if ((do_timings).and.(idx_glbt==-1))       &
         & idx_glbt = psb_get_timer_idx("BLD_AHAT: glob_transp")
    if ((do_timings).and.(idx_phase2==-1))       &
         & idx_phase2 = psb_get_timer_idx("BLD_AHAT: phase2       ")
    if (do_timings) call psb_tic(idx_phase1)
    if (dump_inp) then
      block
        type(psb_ldspmat_type) :: amglob
        write(aname,'(a,i3.3,a)') 'a-bld-inp-',me,'.mtx'
        call a%print(fname=aname,head='Test ')
        call psb_gather(amglob,a,desc_a,info)
        if (me==psb_root_) then
          write(aname,'(a,i0,a)') 'a-bld-inp-g-',amglob%get_nrows(),'.mtx'
          call amglob%print(fname=aname,head='Test ')
        end if
      end block
    end if

    !
    !  Extract diagonal of A
    !
    call a%cp_to(tcoo1)
!!$    call a%triu(atnd,info)
!!$    call atnd%mv_to(tcoo1)
    nr = tcoo1%get_nrows()
    nc = tcoo1%get_ncols()
    nz = tcoo1%get_nzeros()
    diag = a%get_diag(info)
    call psb_realloc(nc,diag,info)
    call psb_halo(diag,desc_a,info)

    ncols = desc_a%get_local_cols()
    call psb_realloc(ncols,ilv,info)
    do i=1, ncols
      ilv(i) = i
    end do
    call desc_a%l2gip(ilv,info,owned=.false.)
    nr = tcoo1%get_nrows()
    nc = tcoo1%get_ncols()
    nz = tcoo1%get_nzeros()
    call tcoo2%allocate(nr,nc,int(1.25*nz))
    k2 = 0
    !
    ! Build the entries of \^A for matching
    !
    minabs = huge(dzero)
    do k = 1, nz
      ii  = tcoo1%ia(k)
      jj  = tcoo1%ja(k)
      !
      ! Run over only the upper triangle, then transpose.
      ! In theory, A was SPD; in practice, since we
      ! are building the hierarchy with Galerkin products P^T A P,
      ! A will be affected by round-off
      !
      if (ilv(ii)<ilv(jj))  then
        aij = tcoo1%val(k)
        aii = diag(ii)
        ajj = diag(jj)
        wii = w(ii)
        wjj = w(jj)
        edgnrm = aii*(wii**2) + ajj*(wjj**2)
        k2 = k2 + 1
        tcoo2%ia(k2) = ii
        tcoo2%ja(k2) = jj
        if (edgnrm > eps) then
          tcoo2%val(k2) = done - (2*done*aij*wii*wjj)/(aii*(wii**2) + ajj*(wjj**2))
        else
          tcoo2%val(k2) = eps
        end if
      end if
      !else
      !  write(0,*) 'build_ahat: index error :',ii,jj,' : boundaries :',nr,nc
      !end if
    end do
    call tcoo2%set_nzeros(k2)

    if (dump) then
      block
        type(psb_ldspmat_type) :: amglob
        call ahat%cp_from(tcoo2)
        call psb_gather(amglob,ahat,desc_a,info)
        if (me==psb_root_) then
          write(aname,'(a,i0,a)') 'ahat-gu-',desc_a%get_global_rows(),'.mtx'
          call amglob%print(fname=aname,head='Test ')
        end if
        write(aname,'(a,i0,a,i0,a)') 'ahat-gu-',desc_a%get_global_rows(),'-p',me,'.mtx'
        call ahat%print(fname=aname,head='Test ',iv=ilv)
        write(aname,'(a,i0,a,i0,a)') 'ainp-ahat-',desc_a%get_global_rows(),'-p',me,'.mtx'
        call a%print(fname=aname,head='Test ',iv=ilv)
        call psb_gather(amglob,a,desc_a,info)
        if (me==psb_root_) then
          write(aname,'(a,i0,a)') 'ainp-ahat-g-',desc_a%get_global_rows(),'.mtx'
          call amglob%print(fname=aname,head='Test ')
        end if
        write(0,*) 'Done build_ahat'
      end block
    end if
    if (do_timings) call psb_toc(idx_phase1)
    if (do_timings) call psb_tic(idx_glbt)

    call psb_glob_transpose(tcoo2,desc_a,info,atrans=tcoo1)
    if (do_timings) call psb_toc(idx_glbt)
    if (do_timings) call psb_tic(idx_phase2)
    if (dump) then
      block
        type(psb_ldspmat_type) :: amglob
        call atnd%cp_from(tcoo1)
        call psb_gather(amglob,atnd,desc_a,info)
        if (me==psb_root_) then
          write(aname,'(a,i0,a)') 'ahat-gl-',desc_a%get_global_rows(),'.mtx'
          call amglob%print(fname=aname,head='Test ')
        end if
        write(aname,'(a,i0,a,i0,a)') 'ahat-gl-',desc_a%get_global_rows(),'-p',me,'.mtx'
        call ahat%print(fname=aname,head='Test ',iv=ilv)
        write(0,*) 'Done build_ahat'
      end block
    end if

    nz  = tcoo1%get_nzeros()
    nz2 = tcoo2%get_nzeros()
    call tcoo2%ensure_size(nz+nz2)
    tcoo2%ia(nz2+1:nz2+nz) = tcoo1%ia(1:nz)
    tcoo2%ja(nz2+1:nz2+nz) = tcoo1%ja(1:nz)
    tcoo2%val(nz2+1:nz2+nz) = tcoo1%val(1:nz)
    call tcoo2%set_nzeros(nz+nz2)
    call tcoo2%fix(info)
    call tcoo1%free()
    nz = tcoo2%get_nzeros()
    if (unroll_logtrans) then
      minabs = huge(dzero)
      do k = 1, nz
        tcoo2%val(k) = abs(tcoo2%val(k))
        minabs = min(minabs,tcoo2%val(k))
      end do
    else
      minabs = minval(abs(tcoo2%val(1:nz)))
    end if
    call psb_min(ictxt,minabs)
    if (minabs <= dzero) then
      if (me == 0) write(0,*) me, 'Min value for log correction is <=zero! '
      minabs = done
    end if
    if (unroll_logtrans) then
      do k = 1, nz
        ! tcoo2%val has already been subject to abs() above.
        tcoo2%val(k) = log(tcoo2%val(k)/(0.999*minabs))
        if (tcoo2%val(k)<0)  write(0,*) me, 'Warning: negative log output!'
      end do
    else
      tcoo2%val(1:nz) = log(abs(tcoo2%val(1:nz))/(0.999*minabs))
      if (any(tcoo2%val(1:nz)<0))  write(0,*) me, 'Warning: negative log output!'
    end if

    call ahat%mv_from(tcoo2)

    if (do_timings) call psb_toc(idx_phase2)

    if (check_size) then
      nzglob = ahat%get_nzeros()
      call psb_sum(ictxt,nzglob)
      if (me==0) write(0,*) 'From build_ahat: global nzeros ',desc_a%get_global_rows(),nzglob
    end if
    if (debug_symmetry) then
      block
        integer(psb_ipk_) :: nz1, nz2
        real(psb_dpk_)    :: mxv
        call ahat%cp_to(tcoo1)
        call psb_glob_transpose(tcoo1,desc_a,info,atrans=tcoo2)
        nz1 = tcoo1%get_nzeros()
        nz2 = tcoo2%get_nzeros()
        call tcoo1%reallocate(nz1+nz2)
        tcoo1%ia(nz1+1:nz1+nz2)  = tcoo2%ia(1:nz2)
        tcoo1%ja(nz1+1:nz1+nz2)  = tcoo2%ja(1:nz2)
        tcoo1%val(nz1+1:nz1+nz2) = -tcoo2%val(1:nz2)
        call tcoo1%set_nzeros(nz1+nz2)
        call tcoo1%set_dupl(psb_dupl_add_)
        call tcoo1%fix(info)
        nz1 = tcoo1%get_nzeros()
        mxv = maxval(abs(tcoo1%val(1:nz1)))
        call psb_max(ictxt,mxv)
        if (me==0) write(0,*) 'Maximum disp from symmetry:',mxv
      end block
    end if

    if (dump_ahat) then
      block
        type(psb_ldspmat_type) :: amglob
        write(aname,'(a,i3.3,a)') 'ahat-hfa-l-',me,'.mtx'
        call ahat%print(fname=aname,head='Test ')
        call psb_gather(amglob,ahat,desc_a,info)
        if (me==psb_root_) then
          write(aname,'(a,i0,a)') 'ahat-hfa-g-',amglob%get_nrows(),'.mtx'
          call amglob%print(fname=aname,head='Test ')
        end if
        write(0,*) 'Done build_ahat'
      end block
    end if

  end subroutine dbuild_ahat

  subroutine psb_dgtranspose(ain,aout,desc_a,info)
    use psb_base_mod
    implicit none
    type(psb_ldspmat_type), intent(in)  :: ain
    type(psb_ldspmat_type), intent(out) :: aout
    type(psb_desc_type)           :: desc_a
    integer(psb_ipk_), intent(out) :: info

    !
    ! BEWARE: This routine works under the assumption
    ! that the same DESC_A works for both A and A^T, which
    ! essentially means that A has a symmetric pattern.
    !
    type(psb_ldspmat_type) :: atmp, ahalo, aglb
    type(psb_ld_coo_sparse_mat) :: tmpcoo
    type(psb_ld_csr_sparse_mat) :: tmpcsr
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: me, np
    integer(psb_lpk_) :: i, j, k, nrow, ncol
    integer(psb_lpk_), allocatable :: ilv(:)
    character(len=80) :: aname
    logical, parameter :: debug=.false., dump=.false., debug_sync=.false.

    ictxt = desc_a%get_context()
    call psb_info(ictxt,me,np)

    nrow = desc_a%get_local_rows()
    ncol = desc_a%get_local_cols()
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*) 'Start gtranspose '
    end if
    call ain%cscnv(tmpcsr,info)

    if (debug) then
      ilv = [(i,i=1,ncol)]
      call desc_a%l2gip(ilv,info,owned=.false.)
      write(aname,'(a,i3.3,a)') 'atmp-preh-',me,'.mtx'
      call ain%print(fname=aname,head='atmp before haloTest ',iv=ilv)
    end if
    if (dump) then
      call ain%cscnv(atmp,info)
      call psb_gather(aglb,atmp,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'aglob-prehalo.mtx'
        call aglb%print(fname=aname,head='Test ')
      end if
    end if

    !call psb_loc_to_glob(tmpcsr%ja,desc_a,info)
    call atmp%mv_from(tmpcsr)

    if (debug) then
      write(aname,'(a,i3.3,a)') 'tmpcsr-',me,'.mtx'
      call atmp%print(fname=aname,head='tmpcsr ',iv=ilv)
      !call psb_set_debug_level(9999)
    end if

    ! FIXME THIS NEEDS REWORKING
    if (debug) write(0,*) me,' Gtranspose into sphalo :',atmp%get_nrows(),atmp%get_ncols()
    call psb_sphalo(atmp,desc_a,ahalo,info,rowscale=.true.)
    if (debug) write(0,*) me,' Gtranspose from sphalo :',ahalo%get_nrows(),ahalo%get_ncols()
    if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=ahalo)

    if (debug) then
      write(aname,'(a,i3.3,a)') 'ahalo-',me,'.mtx'
      call ahalo%print(fname=aname,head='ahalo after haloTest ',iv=ilv)
      write(aname,'(a,i3.3,a)') 'atmp-h-',me,'.mtx'
      call atmp%print(fname=aname,head='atmp after haloTest ',iv=ilv)
    end if

    if (info == psb_success_) call ahalo%free()

    call atmp%cp_to(tmpcoo)
    call tmpcoo%transp()
    !call psb_glob_to_loc(tmpcoo%ia,desc_a,info,iact='I')
    if (debug) write(0,*) 'Before cleanup:',tmpcoo%get_nzeros()

    j = 0
    do k=1, tmpcoo%get_nzeros()
      if ((tmpcoo%ia(k) > 0).and.(tmpcoo%ja(k)>0)) then
        j = j+1
        tmpcoo%ia(j)  = tmpcoo%ia(k)
        tmpcoo%ja(j)  = tmpcoo%ja(k)
        tmpcoo%val(j) = tmpcoo%val(k)
      end if
    end do
    call tmpcoo%set_nzeros(j)

    if (debug) write(0,*) 'After cleanup:',tmpcoo%get_nzeros()

    call ahalo%mv_from(tmpcoo)
    if (dump) then
      call psb_gather(aglb,ahalo,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran-preclip.mtx'
        call aglb%print(fname=aname,head='Test ')
      end if
    end if


    call ahalo%csclip(aout,info,imax=nrow)

    if (debug) write(0,*) 'After clip:',aout%get_nzeros()

    if (debug_sync) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*) 'End gtranspose '
    end if
    !call aout%cscnv(info,type='csr')

    if (dump) then
      write(aname,'(a,i3.3,a)') 'atran-',me,'.mtx'
      call aout%print(fname=aname,head='atrans ',iv=ilv)
      call psb_gather(aglb,aout,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran.mtx'
        call aglb%print(fname=aname,head='Test ')
      end if
    end if

  end subroutine psb_dgtranspose

  subroutine psb_dhtranspose(ain,aout,desc_a,info)
    use psb_base_mod
    implicit none
    type(psb_ldspmat_type), intent(in)  :: ain
    type(psb_ldspmat_type), intent(out) :: aout
    type(psb_desc_type)           :: desc_a
    integer(psb_ipk_), intent(out) :: info

    !
    ! BEWARE: This routine works under the assumption
    ! that the same DESC_A works for both A and A^T, which
    ! essentially means that A has a symmetric pattern.
    !
    type(psb_ldspmat_type) :: atmp, ahalo, aglb
    type(psb_ld_coo_sparse_mat) :: tmpcoo, tmpc1, tmpc2, tmpch
    type(psb_ld_csr_sparse_mat) :: tmpcsr
    integer(psb_ipk_) :: nz1, nz2, nzh, nz
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_) :: me, np
    integer(psb_lpk_) :: i, j, k, nrow, ncol, nlz
    integer(psb_lpk_), allocatable :: ilv(:)
    character(len=80) :: aname
    logical, parameter :: debug=.false., dump=.false., debug_sync=.false.

    ictxt = desc_a%get_context()
    call psb_info(ictxt,me,np)

    nrow = desc_a%get_local_rows()
    ncol = desc_a%get_local_cols()
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*) 'Start htranspose '
    end if
    call ain%cscnv(tmpcsr,info)

    if (debug) then
      ilv = [(i,i=1,ncol)]
      call desc_a%l2gip(ilv,info,owned=.false.)
      write(aname,'(a,i3.3,a)') 'atmp-preh-',me,'.mtx'
      call ain%print(fname=aname,head='atmp before haloTest ',iv=ilv)
    end if
    if (dump) then
      call ain%cscnv(atmp,info)
      call psb_gather(aglb,atmp,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'aglob-prehalo.mtx'
        call aglb%print(fname=aname,head='Test ')
      end if
    end if

    !call psb_loc_to_glob(tmpcsr%ja,desc_a,info)
    call atmp%mv_from(tmpcsr)

    if (debug) then
      write(aname,'(a,i3.3,a)') 'tmpcsr-',me,'.mtx'
      call atmp%print(fname=aname,head='tmpcsr ',iv=ilv)
      !call psb_set_debug_level(9999)
    end if

    ! FIXME THIS NEEDS REWORKING
    if (debug) write(0,*) me,' Htranspose into sphalo :',atmp%get_nrows(),atmp%get_ncols()
    if (.true.) then
      call psb_sphalo(atmp,desc_a,ahalo,info, outfmt='coo  ')
      call atmp%mv_to(tmpc1)
      call ahalo%mv_to(tmpch)
      nz1 = tmpc1%get_nzeros()
      call psb_loc_to_glob(tmpc1%ia(1:nz1),desc_a,info,iact='I')
      call psb_loc_to_glob(tmpc1%ja(1:nz1),desc_a,info,iact='I')
      nzh = tmpch%get_nzeros()
      call psb_loc_to_glob(tmpch%ia(1:nzh),desc_a,info,iact='I')
      call psb_loc_to_glob(tmpch%ja(1:nzh),desc_a,info,iact='I')
      nlz = nz1+nzh
      call tmpcoo%allocate(ncol,ncol,nlz)
      tmpcoo%ia(1:nz1) = tmpc1%ia(1:nz1)
      tmpcoo%ja(1:nz1) = tmpc1%ja(1:nz1)
      tmpcoo%val(1:nz1) = tmpc1%val(1:nz1)
      tmpcoo%ia(nz1+1:nz1+nzh) = tmpch%ia(1:nzh)
      tmpcoo%ja(nz1+1:nz1+nzh) = tmpch%ja(1:nzh)
      tmpcoo%val(nz1+1:nz1+nzh) = tmpch%val(1:nzh)
      call tmpcoo%set_nzeros(nlz)
      call tmpcoo%transp()
      nz = tmpcoo%get_nzeros()
      call psb_glob_to_loc(tmpcoo%ia(1:nz),desc_a,info,iact='I')
      call psb_glob_to_loc(tmpcoo%ja(1:nz),desc_a,info,iact='I')
      if (.true.) then
        call tmpcoo%clean_negidx(info)
      else
        j = 0
        do k=1, tmpcoo%get_nzeros()
          if ((tmpcoo%ia(k) > 0).and.(tmpcoo%ja(k)>0)) then
            j = j+1
            tmpcoo%ia(j)  = tmpcoo%ia(k)
            tmpcoo%ja(j)  = tmpcoo%ja(k)
            tmpcoo%val(j) = tmpcoo%val(k)
          end if
        end do
        call tmpcoo%set_nzeros(j)
      end if
      call ahalo%mv_from(tmpcoo)
      call ahalo%csclip(aout,info,imax=nrow)

    else
      call psb_sphalo(atmp,desc_a,ahalo,info, rowscale=.true.)
      if (debug) write(0,*) me,' Htranspose from sphalo :',ahalo%get_nrows(),ahalo%get_ncols()
      if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=ahalo)

      if (debug) then
        write(aname,'(a,i3.3,a)') 'ahalo-',me,'.mtx'
        call ahalo%print(fname=aname,head='ahalo after haloTest ',iv=ilv)
        write(aname,'(a,i3.3,a)') 'atmp-h-',me,'.mtx'
        call atmp%print(fname=aname,head='atmp after haloTest ',iv=ilv)
      end if

      if (info == psb_success_) call ahalo%free()

      call atmp%cp_to(tmpcoo)
      call tmpcoo%transp()
      if (debug) write(0,*) 'Before cleanup:',tmpcoo%get_nzeros()
      if (.true.) then
        call tmpcoo%clean_negidx(info)
      else

        j = 0
        do k=1, tmpcoo%get_nzeros()
          if ((tmpcoo%ia(k) > 0).and.(tmpcoo%ja(k)>0)) then
            j = j+1
            tmpcoo%ia(j)  = tmpcoo%ia(k)
            tmpcoo%ja(j)  = tmpcoo%ja(k)
            tmpcoo%val(j) = tmpcoo%val(k)
          end if
        end do
        call tmpcoo%set_nzeros(j)
      end if

      if (debug) write(0,*) 'After cleanup:',tmpcoo%get_nzeros()

      call ahalo%mv_from(tmpcoo)
      if (dump) then
        call psb_gather(aglb,ahalo,desc_a,info)
        if (me==psb_root_) then
          write(aname,'(a,i3.3,a)') 'atran-preclip.mtx'
          call aglb%print(fname=aname,head='Test ')
        end if
      end if


      call ahalo%csclip(aout,info,imax=nrow)
    end if

    if (debug) write(0,*) 'After clip:',aout%get_nzeros()

    if (debug_sync) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*) 'End htranspose '
    end if
    !call aout%cscnv(info,type='csr')

    if (dump) then
      write(aname,'(a,i3.3,a)') 'atran-',me,'.mtx'
      call aout%print(fname=aname,head='atrans ',iv=ilv)
      call psb_gather(aglb,aout,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran.mtx'
        call aglb%print(fname=aname,head='Test ')
      end if
    end if

  end subroutine psb_dhtranspose

  subroutine dPMatchBox(nlver,nledge,verlocptr,verlocind,edgelocweight,&
       & verdistance, mate, myrank, numprocs, ictxt,&
       & msgindsent,msgactualsent,msgpercent,&
       & ph0_time, ph1_time, ph2_time, ph1_card, ph2_card,info,display_inp)
    use psb_base_mod
    implicit none
    type(psb_ctxt_type) :: ictxt
    integer(psb_c_ipk_), value :: myrank, numprocs
    integer(psb_c_lpk_), value :: nlver,nledge
    integer(psb_c_lpk_) :: verlocptr(:),verlocind(:), verdistance(:)
    integer(psb_c_lpk_) :: mate(:)
    integer(psb_c_lpk_) :: msgindsent(*),msgactualsent(*)
    real(c_double)      :: ph0_time, ph1_time, ph2_time
    integer(psb_c_lpk_) :: ph1_card(*),ph2_card(*)
    real(c_double)  :: edgelocweight(:)
    real(c_double)      :: msgpercent(*)
    integer(psb_ipk_)   :: info, me, np
    integer(psb_c_mpk_) :: icomm, mrank, mnp
    logical, optional   :: display_inp
    !
    logical, parameter :: debug=.false., debug_out=.false., debug_sync=.false., dump_input=.false.
    logical            :: display_
    integer(psb_lpk_) :: i,k
    integer(psb_ipk_), save :: cnt=1

    call psb_info(ictxt,me,np)
    icomm = psb_get_mpi_comm(ictxt)
    mrank = psb_get_mpi_rank(ictxt,me)
    mnp   = np
    if (present(display_inp)) then
      display_ = display_inp
    else
      display_ = .false.
    end if

    verlocptr(:)   = verlocptr(:)  - 1
    verlocind(:)   = verlocind(:) - 1
    verdistance(:) = verdistance(:) -1

    if (dump_input) then
      block
        integer(psb_ipk_) :: iout=20,info,i,j,k,nr
        character(len=80) :: fname
        write(fname,'(a,i4.4,a,i4.4,a)') 'verlocptr-l',cnt,'-i',me,'.mtx'
        open(iout,file=fname,iostat=info)
        if (info == 0) then
          write(iout,'(a)') '%verlocptr '
          write(iout,*) nlver
          do i=1, nlver+1
            write(iout,*) verlocptr(i)
          end do
          close(iout)
        else
          write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
        end if
        write(fname,'(a,i4.4,a,i4.4,a)') 'edges-l',cnt,'-i',me,'.mtx'
        open(iout,file=fname,iostat=info)
        if (info == 0) then
          write(iout,'(a)') '%verlocind/edgeweight '
          write(iout,*) nlver
          do i=1, nlver
            do j=verlocptr(i)+1,verlocptr(i+1)
              write(iout,*) i-1, verlocind(j),edgelocweight(j)
            end do
          end do
          close(iout)
        else
          write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
        end if
        if (me==0) then
          write(fname,'(a,i4.4,a,i4.4,a)') 'verdistance-l',cnt,'-i',me,'.mtx'
          open(iout,file=fname,iostat=info)
          if (info == 0) then
            write(iout,'(a)') '%verdistance'
            write(iout,*) np
            do i=1, np+1
              write(iout,*) verdistance(i)
            end do
            close(iout)
          else
            write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
          end if
        end if
      end block
      cnt = cnt + 1
    end if

    if (debug.or.display_) then
      do i=0,np-1
        if (me == i) then
          write(6,*) 'Process: ',me,' :  Input into matching: nlver ',nlver,' nledge ',nledge
          write(6,*) 'Process: ',me,' :  VERDISTANCE 0-base : ',verdistance(1:np+1)
          write(6,*) 'Process: ',me,' :  VERLOCPTR   0-base : ',verlocptr(1:nlver+1)
          write(6,*) 'Process: ',me,' :  VERLOCIND   0-base : ',verlocind(1:nledge)
          write(6,*) 'Process: ',me,' :  EDGELOCWEIGHT      : ',edgelocweight(1:nledge)
          write(6,*) 'Process: ',me,' :  Initial MATE       : ',mate(1:nlver)
          flush(6)
        end if
        call psb_barrier(ictxt)
      end do
    end if
    if (debug_sync) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*)' Calling MatchBoxP '
    end if

    call MatchBoxPC(nlver,nledge,verlocptr,verlocind,edgelocweight,&
         & verdistance, mate, mrank, mnp, icomm,&
         & msgindsent,msgactualsent,msgpercent,&
         & ph0_time, ph1_time, ph2_time, ph1_card, ph2_card)
    verlocptr(:)   = verlocptr(:)  + 1
    verlocind(:)   = verlocind(:) + 1
    verdistance(:) = verdistance(:) + 1

    if (debug_sync) then
      call psb_barrier(ictxt)
      if (me == 0) write(0,*)' Done MatchBoxP '
    end if

    if (debug_out) then
      do k=0,np-1
        if (me == k) then
          write(6,*) 'Process: ',me,' :  from Matching (0-base): ',info
          do i=1,nlver
            write(6,*) '(',i,',',mate(i),')'
            !mate(i) = mate(i) +1
            !
          end do
          flush(6)
        end if
        call psb_barrier(ictxt)
      end do
    end if
    where(mate>=0) mate = mate + 1

  end subroutine dPMatchBox

end module dmatchboxp_mod

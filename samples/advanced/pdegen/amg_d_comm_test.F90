!
! Comparison of the PSBLAS communication schemes over an AMG hierarchy.
!
! Companion to test/comm/cg in PSBLAS, which does the same for an unpreconditioned
! CG: here the solve runs with a multilevel preconditioner, so halo exchanges
! happen on every level of the hierarchy and not only on the fine one.
!
! Two measurement modes, both produced by a single invocation so that a batch
! allocation yields the whole dataset:
!
!   uniform      the same scheme on every level, one run per scheme. Answers
!                "which scheme is best overall".
!   sensitivity  baseline everywhere except one level, swept over levels and
!                schemes. Answers "what does scheme S buy me AT level k", which
!                is the quantity a per-level policy is built on and which cannot
!                be recovered from the uniform runs.
!
! The scheme cannot simply be set on desc_a and left at that. Two facts on the
! PSBLAS side make the ordering in run_one mandatory:
!
!   * psb_desc_type%comm_type defaults to psb_comm_isend_irecv_ and is NOT
!     carried over by clone/transfer, so the coarse descriptors built inside
!     prec%build come out on the baseline scheme whatever desc_a says;
!   * the communication handle is latched lazily on the vector at its first
!     exchange (psi_dswapdata_vect) and never re-read afterwards, so handles
!     left over from a previous repetition must be freed explicitly.
!
! Run with:
!   mpirun -np <P> ./amg_d_comm_test [idim] [nrep] [nlev] [itmax]
!                                    [--remap] [--csv=<path>] [--mode=all|uniform|sensitivity]
!
program amg_d_comm_test
  use psb_base_mod
  use psb_util_mod
  use psb_linsolve_mod
  use amg_prec_mod
  use amg_d_genpde_mod
  use amg_d_pde3d_poisson_mod
  use psb_comm_factory_mod, only: psb_comm_free
  use psb_comm_schemes_mod, only: psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
       & psb_comm_persistent_ineighbor_alltoallv_, psb_comm_rma_pull_, psb_comm_rma_push_

  implicit none

  integer(psb_ipk_), parameter :: n_schemes = 5, max_lev = 32

  type(psb_ctxt_type)   :: ctxt
  type(psb_dspmat_type) :: a
  type(psb_desc_type)   :: desc_a
  type(psb_d_vect_type) :: b, x
  type(amg_dprec_type)  :: prec

  integer(psb_ipk_) :: info, my_rank, np
  integer(psb_ipk_) :: idim, itmax, itrace, istop
  integer(psb_ipk_) :: s_idx, rep, nrep, nlev, lev, nlev_built, csv_unit
  integer(psb_ipk_) :: scheme_type(n_schemes)
  integer(psb_ipk_) :: scheme_of_level(max_lev)
  character(len=25) :: scheme_name(n_schemes)
  real(psb_dpk_)    :: eps
  character(len=5)  :: afmt
  character(len=256):: csv_file, lev_file
  character(len=16) :: run_mode
  logical           :: do_remap, ok_all, want_csv, header_done
  ! per-level exchange accounting, filled after every solve
  integer(psb_ipk_) :: lev_unit
  integer(psb_ipk_) :: lev_rows(max_lev), lev_halo(max_lev)

  info   = psb_success_
  afmt   = 'CSR'
  idim   = 40
  nrep   = 5
  nlev   = 3
  itmax  = 1000
  itrace = -1
  istop  = 2
  eps    = 1.0e-6_psb_dpk_
  csv_file = ''
  run_mode = 'all'
  do_remap = .false.
  csv_unit = 77
  lev_unit = 78
  header_done = .false.
  lev_rows(:)  = 0
  lev_halo(:)  = 0

  scheme_type = (/ psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
       & psb_comm_persistent_ineighbor_alltoallv_, psb_comm_rma_pull_, psb_comm_rma_push_ /)
  scheme_name(1) = 'isend_irecv'
  scheme_name(2) = 'ineighbor_alltoallv'
  scheme_name(3) = 'persistent_ineighbor_a2av'
  scheme_name(4) = 'rma_pull'
  scheme_name(5) = 'rma_push'

  call read_int_arg(1, idim,  40)
  call read_int_arg(2, nrep,   5)
  call read_int_arg(3, nlev,   3)
  call read_int_arg(4, itmax,1000)
  call parse_flags(do_remap, csv_file, run_mode)
  want_csv = (len_trim(csv_file) > 0)

  call psb_init(ctxt)
  call psb_info(ctxt, my_rank, np)
  call amg_set_do_remap(do_remap)

  if (my_rank == psb_root_) then
    write(psb_out_unit,*) 'Welcome to PSBLAS version: ', psb_version_string_
    write(psb_out_unit,'("AMG communication-scheme comparison")')
    write(psb_out_unit,'("Grid dimensions      : ",i0," x ",i0," x ",i0)') idim,idim,idim
    write(psb_out_unit,'("Number of processors : ",i0)') np
    write(psb_out_unit,'("Preconditioner       : ML, V-cycle, JACOBI, max levels ",i0)') nlev
    write(psb_out_unit,'("Iterative method     : CG, itmax ",i0,", eps ",es9.2)') itmax, eps
    write(psb_out_unit,'("Repetitions          : ",i0)') nrep
    write(psb_out_unit,'("Remap active         : ",l1)') amg_get_do_remap()
    write(psb_out_unit,'("Mode                 : ",a)') trim(run_mode)
    if (want_csv) write(psb_out_unit,'("CSV output           : ",a)') trim(csv_file)
    write(psb_out_unit,'(" ")')
  end if

  if (want_csv .and. (my_rank == psb_root_)) then
    open(unit=csv_unit,file=trim(csv_file),status='replace',action='write',iostat=info)
    if (info /= 0) then
      write(psb_err_unit,'("Cannot open CSV file ",a)') trim(csv_file)
      goto 9999
    end if
    write(csv_unit,'(a)') 'mode,scheme,target_level,nranks,idim,nlev_built,rep,'// &
         & 'prec_init_s,prec_bld_s,comm_set_s,krylov_s,total_s,iters,final_err,remap'
    ! Second file, one row per (run, level): this is the per-level exchange
    ! breakdown the policy is meant to be built on.
    lev_file = trim(csv_file)//'.levels.csv'
    open(unit=lev_unit,file=trim(lev_file),status='replace',action='write',iostat=info)
    if (info /= 0) then
      write(psb_err_unit,'("Cannot open CSV file ",a)') trim(lev_file)
      goto 9999
    end if
    write(lev_unit,'(a)') 'mode,scheme,target_level,nranks,idim,rep,level,'// &
         & 'loc_rows,halo_width'
    header_done = .true.
  end if

  call psb_barrier(ctxt)
  call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
       & a1_poisson,a2_poisson,a3_poisson,&
       & b1_poisson,b2_poisson,b3_poisson,c_poisson,g_poisson,info)
  if (info /= psb_success_) goto 9999

  ok_all = .true.
  nlev_built = 0

  ! ---------------------------------------------------------------- uniform
  if ((trim(run_mode) == 'all').or.(trim(run_mode) == 'uniform')) then
    if (my_rank == psb_root_) then
      write(psb_out_unit,'(104("="))')
      write(psb_out_unit,'("UNIFORM: same scheme on every level")')
      write(psb_out_unit,'(104("="))')
    end if
    do s_idx = 1, n_schemes
      scheme_of_level(:) = scheme_type(s_idx)
      do rep = 1, nrep
        call run_one('uniform', trim(scheme_name(s_idx)), izero, &
             & scheme_of_level, (rep == 1), info)
        if (info /= psb_success_) goto 9999
      end do
    end do
  end if

  ! ------------------------------------------------------------ sensitivity
  !
  ! Baseline everywhere, one level at a time moved onto another scheme. The
  ! difference against the all-baseline uniform run is the marginal value of
  ! that scheme at that level, which is exactly what a per-level policy needs.
  !
  if ((trim(run_mode) == 'all').or.(trim(run_mode) == 'sensitivity')) then
    if (nlev_built <= 0) then
      ! Need one build to learn how many levels the hierarchy actually has.
      scheme_of_level(:) = scheme_type(1)
      call run_one('probe', trim(scheme_name(1)), izero, scheme_of_level, .false., info)
      if (info /= psb_success_) goto 9999
    end if
    if (my_rank == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'(104("="))')
      write(psb_out_unit,'("SENSITIVITY: baseline everywhere except one level")')
      write(psb_out_unit,'(104("="))')
    end if
    do lev = 1, nlev_built
      do s_idx = 2, n_schemes
        scheme_of_level(:)   = scheme_type(1)
        scheme_of_level(lev) = scheme_type(s_idx)
        do rep = 1, nrep
          call run_one('sensitivity', trim(scheme_name(s_idx)), lev, &
               & scheme_of_level, .false., info)
          if (info /= psb_success_) goto 9999
        end do
      end do
    end do
  end if

  if (my_rank == psb_root_) then
    write(psb_out_unit,'(" ")')
    if (ok_all) then
      write(psb_out_unit,'("SCHEME PROPAGATION: OK on every level and every configuration")')
    else
      write(psb_out_unit,'("SCHEME PROPAGATION: FAILED, see the messages above")')
    end if
    if (want_csv) then
      close(csv_unit)
      close(lev_unit)
      write(psb_out_unit,'("CSV written to ",a)') trim(csv_file)
    end if
    write(psb_out_unit,'(104("="))')
  end if

  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call psb_cdfree(desc_a,info)
  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)
  stop 1

contains

  !
  ! One full init/build/solve with the given per-level scheme assignment.
  ! Timings are max-reduced across ranks before being recorded.
  !
  subroutine run_one(mode, sname, target_lev, sch_of_lev, verbose, info)
    character(len=*), intent(in)   :: mode, sname
    integer(psb_ipk_), intent(in)  :: target_lev
    integer(psb_ipk_), intent(in)  :: sch_of_lev(:)
    logical, intent(in)            :: verbose
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_) :: t_start, t_init, t_bld, t_comm, t_kry, err
    integer(psb_ipk_) :: iter, lv

    info = psb_success_

    ! Fine level: must be set before the build, so that desc_a and every vector
    ! lazily initialised during the solve pick this scheme up.
    call desc_a%set_comm_scheme(sch_of_lev(1), info)
    if (info /= psb_success_) return
    if (allocated(x%v%comm_handle)) call psb_comm_free(x%v%comm_handle, info)
    if (info /= psb_success_) return
    if (allocated(b%v%comm_handle)) call psb_comm_free(b%v%comm_handle, info)
    if (info /= psb_success_) return

    call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
    if (info /= psb_success_) return
    call psb_barrier(ctxt)

    t_start = psb_wtime()
    call prec%init(ctxt,'ML',info)
    if (info /= psb_success_) return
    ! Spelled out rather than left to the defaults, so that the configuration
    ! the timings refer to is readable here and not in amg_dprecinit.
    call prec%set('max_levs',      nlev,     info)
    if (info == psb_success_) call prec%set('ml_cycle','VCYCLE',info)
    if (info == psb_success_) call prec%set('smoother_type','JACOBI',info)
    if (info == psb_success_) call prec%set('smoother_sweeps',ione,info)
    if (info /= psb_success_) return
    t_init = psb_wtime() - t_start
    call psb_amx(ctxt,t_init)

    t_start = psb_wtime()
    call prec%build(a,desc_a,info)
    if (info /= psb_success_) return
    t_bld = psb_wtime() - t_start
    call psb_amx(ctxt,t_bld)

    nlev_built = prec%get_nlevs()

    ! Hierarchy: the coarse descriptors only exist now, and they were born on
    ! the default scheme. This must happen before the first apply.
    t_start = psb_wtime()
    call set_hierarchy_scheme(prec, sch_of_lev, info)
    if (info /= psb_success_) return
    t_comm = psb_wtime() - t_start
    call psb_amx(ctxt,t_comm)

    call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
    if (info /= psb_success_) return
    call psb_barrier(ctxt)

    t_start = psb_wtime()
    call psb_krylov('CG',a,prec,b,x,eps,desc_a,info,&
         & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istop)
    t_kry = psb_wtime() - t_start
    call psb_amx(ctxt,t_kry)
    if (info /= psb_success_) return

    call check_hierarchy_scheme(prec, sch_of_lev, sname, my_rank, ok_all)

    ! Per-level structure. The exchange time per level would need counters
    ! inside psi_swapdata: adding any component to psb_desc_type or to
    ! psb_comm_handle_type makes gfortran 13.3 die with an ICE in
    ! amg_?_matchboxp_mod (gfc_get_function_type), so that is left out and the
    ! per-level question is answered by the sensitivity mode instead.
    do lv = 1, nlev_built
      lev_rows(lv) = 0
      lev_halo(lv) = 0
      if (prec%precv(lv)%desc_ac%is_asb()) then
        lev_rows(lv) = prec%precv(lv)%desc_ac%get_local_rows()
        lev_halo(lv) = prec%precv(lv)%desc_ac%get_local_cols() - lev_rows(lv)
      end if
    end do

    if (my_rank == psb_root_) then
      if (verbose) call report_hierarchy(prec)
      if (want_csv .and. header_done .and. (trim(mode) /= 'probe')) then
        do lv = 1, nlev_built
          write(lev_unit,'(a,",",a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0)') &
               & trim(mode), trim(sname), target_lev, np, idim, rep, lv, &
               & lev_rows(lv), lev_halo(lv)
        end do
      end if
      if (trim(mode) /= 'probe') then
        write(psb_out_unit,'(a14,1x,a26,1x,"lev ",i2,1x,"it ",i5,1x,"err ",es12.5,&
             &1x,"bld ",es11.4,1x,"solve ",es11.4)') &
             & trim(mode), trim(sname), target_lev, iter, err, t_bld, t_kry
      end if
      if (want_csv .and. header_done .and. (trim(mode) /= 'probe')) then
        write(csv_unit,'(a,",",a,",",i0,",",i0,",",i0,",",i0,",",i0,5(",",es16.9),",",i0,",",es16.9,",",l1)') &
             & trim(mode), trim(sname), target_lev, np, idim, nlev_built, rep, &
             & t_init, t_bld, t_comm, t_kry, t_init+t_bld+t_comm+t_kry, &
             & iter, err, do_remap
      end if
    end if

    call prec%free(info)
  end subroutine run_one

  !
  ! Write the per-level schemes onto every descriptor of the hierarchy and drop
  ! the handles already latched on the per-level work vectors.
  !
  ! base_desc is a pointer (level 1 aliases desc_a); desc_ac is the coarse
  ! descriptor produced by the aggregation during prec%build. Both must be set:
  ! the smoother exchanges go through base_desc, the coarse-grid ones through
  ! desc_ac.
  !
  subroutine set_hierarchy_scheme(prec, sch_of_lev, info)
    type(amg_dprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in)       :: sch_of_lev(:)
    integer(psb_ipk_), intent(out)      :: info
    integer(psb_ipk_) :: lv, iv, sch

    info = psb_success_
    if (.not.allocated(prec%precv)) return

    do lv = 1, size(prec%precv)
      sch = sch_of_lev(min(lv,size(sch_of_lev)))
      if (associated(prec%precv(lv)%base_desc)) then
        call prec%precv(lv)%base_desc%set_comm_scheme(sch, info)
        if (info /= psb_success_) return
      end if
      call prec%precv(lv)%desc_ac%set_comm_scheme(sch, info)
      if (info /= psb_success_) return

      if (allocated(prec%precv(lv)%wrk)) then
        call free_vect_handle(prec%precv(lv)%wrk%vtx,  info)
        if (info == psb_success_) call free_vect_handle(prec%precv(lv)%wrk%vty,  info)
        if (info == psb_success_) call free_vect_handle(prec%precv(lv)%wrk%vx2l, info)
        if (info == psb_success_) call free_vect_handle(prec%precv(lv)%wrk%vy2l, info)
        if (info /= psb_success_) return
        if (allocated(prec%precv(lv)%wrk%wv)) then
          do iv = 1, size(prec%precv(lv)%wrk%wv)
            call free_vect_handle(prec%precv(lv)%wrk%wv(iv), info)
            if (info /= psb_success_) return
          end do
        end if
      end if
    end do
  end subroutine set_hierarchy_scheme

  subroutine free_vect_handle(v, info)
    type(psb_d_vect_type), intent(inout) :: v
    integer(psb_ipk_), intent(out)       :: info
    info = psb_success_
    if (.not.allocated(v%v)) return
    if (allocated(v%v%comm_handle)) call psb_comm_free(v%v%comm_handle, info)
  end subroutine free_vect_handle

  !
  ! After a solve, every work vector that took part in an exchange must carry a
  ! handle of the requested type. A level still holding the baseline while the
  ! others moved is the signature of an incomplete propagation.
  !
  subroutine check_hierarchy_scheme(prec, sch_of_lev, sname, my_rank, ok_all)
    type(amg_dprec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in)    :: sch_of_lev(:), my_rank
    character(len=*), intent(in)     :: sname
    logical, intent(inout)           :: ok_all
    integer(psb_ipk_) :: lv, got, want

    if (.not.allocated(prec%precv)) return
    do lv = 1, size(prec%precv)
      if (.not.allocated(prec%precv(lv)%wrk)) cycle
      if (.not.allocated(prec%precv(lv)%wrk%vtx%v)) cycle
      if (.not.allocated(prec%precv(lv)%wrk%vtx%v%comm_handle)) cycle
      want = sch_of_lev(min(lv,size(sch_of_lev)))
      got  = prec%precv(lv)%wrk%vtx%v%comm_handle%comm_type
      if (got /= want) then
        ok_all = .false.
        write(psb_err_unit,'("  [FAIL] ",a,": rank ",i0," level ",i0,": expected ",i0," got ",i0)') &
             & sname, my_rank, lv, want, got
      end if
    end do
  end subroutine check_hierarchy_scheme

  !
  ! Per-level size and halo width. The halo is the gap between local columns and
  ! local rows, i.e. how many remote entries this rank receives: it is the
  ! quantity whose ratio to the local rows degrades going coarse, which is the
  ! regime the schemes are expected to behave differently in.
  !
  subroutine report_hierarchy(prec)
    type(amg_dprec_type), intent(in) :: prec
    integer(psb_ipk_) :: lv, nr, nc

    if (.not.allocated(prec%precv)) return
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("  level      loc rows      loc cols     halo width   halo/rows")')
    do lv = 1, size(prec%precv)
      if (.not.prec%precv(lv)%desc_ac%is_asb()) cycle
      nr = prec%precv(lv)%desc_ac%get_local_rows()
      nc = prec%precv(lv)%desc_ac%get_local_cols()
      write(psb_out_unit,'(i7,i14,i14,i15,f12.3)') lv, nr, nc, nc-nr, &
           & real(nc-nr,psb_dpk_)/real(max(nr,1),psb_dpk_)
    end do
    write(psb_out_unit,'(" ")')
  end subroutine report_hierarchy

  subroutine parse_flags(do_remap, csv_file, run_mode)
    logical, intent(inout)            :: do_remap
    character(len=*), intent(inout)   :: csv_file, run_mode
    character(len=256) :: larg
    integer(psb_ipk_)  :: k

    do k = 1, command_argument_count()
      call get_command_argument(k, larg)
      if (index(larg,'--remap') == 1) then
        do_remap = .true.
      else if (index(larg,'--csv=') == 1) then
        csv_file = trim(larg(7:))
      else if (index(larg,'--mode=') == 1) then
        run_mode = trim(larg(8:))
      end if
    end do
  end subroutine parse_flags

  subroutine read_int_arg(pos, val, fallback)
    integer(psb_ipk_), intent(in)    :: pos, fallback
    integer(psb_ipk_), intent(inout) :: val
    character(len=256) :: larg
    integer(psb_ipk_)  :: ios

    call get_command_argument(pos, larg)
    if (len_trim(larg) <= 0) return
    if (index(larg,'--') == 1) return
    read(larg,*,iostat=ios) val
    if ((ios /= 0).or.(val <= 0)) val = fallback
  end subroutine read_int_arg

end program amg_d_comm_test

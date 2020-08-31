!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.0)
!    
!    (C) Copyright 2008,2009,2010,2012,2013
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
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!
!  Current version of this file contributed by:
!        Ambra Abdullahi Hassan 
!
!  
subroutine d_mumps_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

  use psb_base_mod
  use amg_d_mumps_solver
  Implicit None

  ! Arguments
  type(psb_dspmat_type)                               :: c 
  type(psb_dspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(inout)                  :: desc_a 
  class(amg_d_mumps_solver_type), intent(inout)       :: sv
  integer(psb_ipk_), intent(out)                      :: info
  type(psb_dspmat_type), intent(in), target, optional :: b
  class(psb_d_base_sparse_mat), intent(in), optional  :: amold
  class(psb_d_base_vect_type), intent(in), optional   :: vmold
  class(psb_i_base_vect_type), intent(in), optional   :: imold
  ! Local variables
  type(psb_dspmat_type)      :: atmp
  type(psb_d_coo_sparse_mat), target :: acoo
#if defined(IPK4) && defined(LPK8)
  integer(psb_lpk_), allocatable :: gia(:), gja(:)
#endif
  integer(psb_ipk_)  :: n_row,n_col, nrow_a, nza, npr, npc
  integer(psb_lpk_)  :: nglob, nglobrec, nzt
  integer(psb_ipk_)  :: ifrst, ibcheck
  integer(psb_ipk_)  :: ictxt, ictxt1, icomm, np, iam, me, i, err_act, debug_unit, debug_level
  character(len=20)  :: name='d_mumps_solver_bld', ch_err

#if defined(HAVE_MUMPS_) 

  info=psb_success_

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, iam, np)      
  if (sv%ipar(1) == amg_local_solver_ ) then
    call psb_init(ictxt1,np=1,basectxt=ictxt,ids=(/iam/))
    icomm = psb_get_mpi_comm(ictxt1)
    allocate(sv%local_ictxt,stat=info)
    sv%local_ictxt = ictxt1
    !write(*,*)iam,'mumps_bld: local +++++>',icomm,sv%local_ictxt
    call psb_info(ictxt1, me, np)
    npr  = np
  else if (sv%ipar(1) == amg_global_solver_ ) then
    icomm = psb_get_mpi_comm(ictxt)
    !write(*,*)iam,'mumps_bld: global +++++>',icomm,ictxt
    call psb_info(ictxt, iam, np)
    me = iam 
    npr  = np
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,& 
         & a_err='Invalid local/global solver in MUMPS')
    goto 9999
  end if
  npc  = 1
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'
  !    if (allocated(sv%id)) then 
  !      call sv%free(info)

  !     deallocate(sv%id)
  !   end if
  if(.not.allocated(sv%id)) then
    allocate(sv%id,stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name,a_err='amg_dmumps_default')
      goto 9999
    end if
  end if


  sv%id%comm    =  icomm
  sv%id%job     = -1
  sv%id%par     =  1
  if (sv%ipar(3) == 2) then
    sv%id%sym = 2
  else
    sv%id%sym = 0
  end if

  call dmumps(sv%id)   
  !WARNING: CALLING dmumps WITH JOB=-1 DESTROY THE SETTING OF DEFAULT:TO FIX
  if (allocated(sv%icntl)) then 
    do i=1,amg_mumps_icntl_size
      if (allocated(sv%icntl(i)%item)) then
        !write(0,*) 'MUMPS_BLD: setting entry ',i,' to ', sv%icntl(i)%item
        sv%id%icntl(i) = sv%icntl(i)%item
      end if
    end do
  end if
  if (allocated(sv%rcntl)) then 
    do i=1,amg_mumps_rcntl_size
      if (allocated(sv%rcntl(i)%item)) sv%id%cntl(i) = sv%rcntl(i)%item
    end do
  end if
  sv%id%icntl(3)=sv%ipar(2)

  nglob  = desc_a%get_global_rows()
  if (sv%ipar(1) == amg_local_solver_ ) then
    nglobrec=desc_a%get_local_rows()
    if (sv%ipar(3) == 2) then
      ! Always pass the upper triangle to MUMPS
      call a%triu(c,info,jmax=a%get_nrows())
      call c%set_symmetric()
    else
      call a%csclip(c,info,jmax=a%get_nrows())
    end if
    call c%cp_to(acoo)
    nglob = c%get_nrows()
    if (nglobrec /= nglob) then
      write(*,*)'WARNING: MUMPS solver does not allow overlap in AS yet. '
      write(*,*)'A zero-overlap is used instead'
    end if
  else
    call a%cp_to(acoo)
  end if
  nza = acoo%get_nzeros()

  ! switch to global numbering
  if (sv%ipar(1) == amg_global_solver_ ) then
#if defined(IPK4) && defined(LPK8)
    !
    ! Strategy here is as follows: because a call to MUMPS
    ! as a gobal solver is mostly done  at the coarsest level,
    ! even if we start from a problem requiring 8 bytes, chances
    ! are that the global size will be suitable for 4 bytes
    ! anyway, so we hope for the best, and throw an error
    ! if something goes wrong.
    !
    if (nglob > huge(1_psb_ipk_)) then
      write(0,*) iam,' ',trim(name),': Error: overflow of local indices '
      info=psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999
    end if

    gia = acoo%ia(1:nza)
    gja = acoo%ja(1:nza)
    call psb_loc_to_glob(gia(1:nza), desc_a, info, iact='I')
    call psb_loc_to_glob(gja(1:nza), desc_a, info, iact='I')
    acoo%ia(1:nza) = gia(1:nza)
    acoo%ja(1:nza) = gja(1:nza)
#else
    !
    ! Here global and local numbers have the same size, so this must work.
    !
    call psb_loc_to_glob(acoo%ja(1:nza), desc_a, info, iact='I')
    call psb_loc_to_glob(acoo%ia(1:nza), desc_a, info, iact='I')
#endif
    if (sv%ipar(3) == 2 ) then
      !  Always pass the upper triangle to MUMPS
      block
        integer(psb_ipk_) :: j,nz
        nz = 0
        do j=1,nza
          if (acoo%ja(j) >= acoo%ia(j)) then
            nz = nz + 1
            acoo%ia(nz)  = acoo%ia(j)
            acoo%ja(nz)  = acoo%ja(j)
            acoo%val(nz) = acoo%val(j)
          end if
        end do
        call acoo%set_nzeros(nz)
        call acoo%set_triangle()
        call acoo%set_upper()
        call acoo%set_symmetric()
      end block
    end if
  end if
  sv%id%irn_loc   => acoo%ia
  sv%id%jcn_loc   => acoo%ja
  sv%id%a_loc     => acoo%val
  sv%id%icntl(18) = 3
  sv%id%n      = nglob
  ! there should be a better way for this
  sv%id%nnz_loc = acoo%get_nzeros()
  sv%id%nnz     = acoo%get_nzeros()
  sv%id%job    = 4
  if (sv%ipar(1) == amg_global_solver_ ) then
    call psb_sum(ictxt,sv%id%nnz)
  end if
  !call psb_barrier(ictxt)
  write(*,*)iam, ' calling mumps N,nz,nz_loc',sv%id%n,sv%id%nnz,sv%id%nnz_loc
  call dmumps(sv%id)
  !call psb_barrier(ictxt)
  info = sv%id%infog(1)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='amg_dmumps_fact '
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  nullify(sv%id%irn)
  nullify(sv%id%jcn)
  nullify(sv%id%a)

  call acoo%free()
  sv%built=.true.

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) iam,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
#else
  write(psb_err_unit,*) "MUMPS Not Configured, fix make.inc and recompile "
#endif
end subroutine d_mumps_solver_bld


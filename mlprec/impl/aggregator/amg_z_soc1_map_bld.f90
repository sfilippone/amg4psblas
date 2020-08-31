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
!
! File: amg_z_soc1_map__bld.f90
!
! Subroutine: amg_z_soc1_map_bld
! Version:    complex
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
!    a       -  type(psb_zspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_zprec_type), input/output.
!               The preconditioner data structure; upon exit it contains 
!               the multilevel hierarchy of prolongators, restrictors
!               and coarse matrices.
!    info    -  integer, output.
!               Error code.              
!
!
!
subroutine amg_z_soc1_map_bld(iorder,theta,clean_zeros,a,desc_a,nlaggr,ilaggr,info)

  use psb_base_mod
  use amg_base_prec_type
  use amg_z_inner_mod

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)     :: iorder
  logical, intent(in)               :: clean_zeros
  type(psb_zspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)    :: desc_a
  real(psb_dpk_), intent(in)         :: theta
  integer(psb_lpk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
  integer(psb_ipk_), intent(out)               :: info

  ! Local variables
  integer(psb_ipk_), allocatable :: ils(:), neigh(:), irow(:), icol(:),&
       & ideg(:), idxs(:)
  integer(psb_lpk_), allocatable :: tmpaggr(:)
  complex(psb_dpk_), allocatable  :: val(:), diag(:)
  integer(psb_ipk_) :: icnt,nlp,k,n,ia,isz,nr, nc, naggr,i,j,m, nz, ilg, ii, ip
  type(psb_z_csr_sparse_mat) :: acsr
  real(psb_dpk_)  :: cpling, tcl
  logical :: disjoint
  integer(psb_ipk_) :: debug_level, debug_unit,err_act
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: nrow, ncol, n_ne
  integer(psb_lpk_) :: nrglob
  character(len=20)  :: name, ch_err

  info=psb_success_
  name = 'amg_soc1_map_bld'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  !
  ictxt=desc_a%get_context()
  call psb_info(ictxt,me,np)
  nrow   = desc_a%get_local_rows()
  ncol   = desc_a%get_local_cols()
  nrglob = desc_a%get_global_rows()

  nr = a%get_nrows()
  nc = a%get_ncols()
  allocate(ilaggr(nr),neigh(nr),ideg(nr),idxs(nr),&
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

  call a%cp_to(acsr)
  if (clean_zeros) call acsr%clean_zeros(info)
  if (iorder == amg_aggr_ord_nat_) then 
    do i=1, nr
      ilaggr(i) = -(nr+1)
      idxs(i)   = i 
    end do
  else 
    do i=1, nr
      ilaggr(i) = -(nr+1)
      ideg(i)   = acsr%irp(i+1) - acsr%irp(i)
    end do
    call psb_msort(ideg,ix=idxs,dir=psb_sort_down_)
  end if


  !
  ! Phase one: Start with disjoint groups.
  ! 
  naggr = 0
  icnt = 0
  step1: do ii=1, nr
    i = idxs(ii)
    if ((i<1).or.(i>nr)) then
      info=psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    
    if (ilaggr(i) == -(nr+1)) then
      nz         = (acsr%irp(i+1)-acsr%irp(i))
      if ((nz<0).or.(nz>size(icol))) then
        info=psb_err_internal_error_
        call psb_errpush(info,name)
        goto 9999
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
        icnt      = icnt + 1 
        naggr     = naggr + 1
        do k=1, ip
          ilaggr(icol(k)) = naggr
        end do
        ilaggr(i) = naggr
      end if
    endif
  enddo step1
  
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
        icnt      = icnt + 1 
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

  if (naggr > ncol) then 
    !write(0,*) name,'Error : naggr > ncol',naggr,ncol
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
  call psb_sum(ictxt,nlaggr(1:np))

  call acsr%free()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_z_soc1_map_bld


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
! File: amg_zaggrmat_minnrg_bld.F90
!
! Subroutine: amg_zaggrmat_minnrg_bld
! Version:    complex
!
!  This routine builds a coarse-level matrix A_C from a fine-level matrix A
!  by using the Galerkin approach, i.e.
!
!                               A_C = P_C^T A P_C,
!
!  where P_C is a prolongator from the coarse level to the fine one.
! 
!  The prolongator P_C is built according to a smoothed aggregation algorithm,
!  i.e. it is obtained by applying a damped Jacobi smoother to the piecewise
!  constant interpolation operator P corresponding to the fine-to-coarse level 
!  mapping built by the amg_aggrmap_bld subroutine:
!
!                            P_C = (I - omega*D^(-1)A) * P,
!
!  where D is the diagonal matrix with main diagonal equal to the main diagonal
!  of A, and omega is a suitable smoothing parameter. An estimate of the spectral
!  radius of D^(-1)A, to be used in the computation of omega, is provided, 
!  according to the value of p%parms%aggr_omega_alg, specified by the user
!  through amg_zprecinit and amg_zprecset.
!  4. Minimum energy aggregation:
!    M. Sala, R. Tuminaro: A new Petrov-Galerkin smoothed aggregation preconditioner
!    for nonsymmetric linear systems, SIAM J. Sci. Comput., 31(1):143-166 (2008)
!     
!  On output from this routine the entries of AC, op_prol, op_restr
!  are still in "global numbering" mode; this is fixed in the calling routine
!  aggregator%mat_bld.
!
!
! Arguments:
!    a          -  type(psb_zspmat_type), input.     
!                  The sparse matrix structure containing the local part of
!                  the fine-level matrix.
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the fine-level matrix.
!    p          -  type(amg_z_onelev_type), input/output.
!                  The 'one-level' data structure that will contain the local
!                  part of the matrix to be built as well as the information
!                  concerning the prolongator and its transpose.
!    parms      -   type(amg_dml_parms), input
!                  Parameters controlling the choice of algorithm
!    ac         -  type(psb_zspmat_type), output
!                  The coarse matrix on output 
!                  
!    ilaggr     -  integer, dimension(:), input
!                  The mapping between the row indices of the coarse-level
!                  matrix and the row indices of the fine-level matrix.
!                  ilaggr(i)=j means that node i in the adjacency graph
!                  of the fine-level matrix is mapped onto node j in the
!                  adjacency graph of the coarse-level matrix. Note that the indices
!                  are assumed to be shifted so as to make sure the ranges on
!                  the various processes do not   overlap.
!    nlaggr     -  integer, dimension(:) input
!                  nlaggr(i) contains the aggregates held by process i.
!    op_prol    -  type(psb_zspmat_type), input/output
!                  The tentative prolongator on input, the computed prolongator on output
!               
!    op_restr    -  type(psb_zspmat_type), output
!                  The restrictor operator; in this particular case, it is different
!                  from the transpose of the prolongator. 
!               
!    info       -  integer, output.
!                  Error code.
!
!
subroutine amg_zaggrmat_minnrg_bld(a,desc_a,ilaggr,nlaggr,parms,&
     & ac,desc_ac,op_prol,op_restr,info)
  use psb_base_mod
  use amg_base_prec_type
  use amg_z_inner_mod, amg_protect_name => amg_zaggrmat_minnrg_bld

  implicit none

  ! Arguments
  type(psb_zspmat_type), intent(in)           :: a
  type(psb_desc_type), intent(inout)            :: desc_a
  integer(psb_lpk_), intent(inout)              :: ilaggr(:), nlaggr(:)
  type(amg_dml_parms), intent(inout)         :: parms 
  type(psb_lzspmat_type), intent(inout)       :: op_prol
  type(psb_lzspmat_type), intent(out)         :: ac,op_restr
  type(psb_desc_type), intent(inout)            :: desc_ac
  integer(psb_ipk_), intent(out)                :: info

  ! Local variables
  integer(psb_lpk_) :: nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
       & naggr, nzl,naggrm1,naggrp1, i, j, k, jd, icolF, nrt
  type(psb_ctxt_type)       :: ctxt
  integer(psb_ipk_)         :: np, me, icomm
  character(len=20)         :: name
  type(psb_lzspmat_type)      :: la, af, ptilde, rtilde, atran, atp, atdatp
  type(psb_lzspmat_type)      :: am3,am4, ap, adap,atmp,rada, ra, atmp2, dap, dadap, da
  type(psb_lzspmat_type)      :: dat, datp, datdatp, atmp3, tmp_prol
  type(psb_lz_coo_sparse_mat) :: tmpcoo
  type(psb_lz_csr_sparse_mat) :: acsr1, acsr2, acsr3, acsr, acsrf
  type(psb_lz_csc_sparse_mat) :: csc_dap, csc_dadap, csc_datp, csc_datdatp, acsc
  complex(psb_dpk_), allocatable :: adiag(:), adinv(:)
  complex(psb_dpk_), allocatable :: omf(:), omp(:), omi(:), oden(:)
  logical                    :: filter_mat
  integer(psb_ipk_)          :: ierr(5)
  integer(psb_ipk_)          :: debug_level, debug_unit, err_act
  integer(psb_ipk_), parameter :: ncmax=16
  real(psb_dpk_)              :: anorm, theta
  complex(psb_dpk_)            :: tmp, alpha, beta, ommx

  name='amg_aggrmat_minnrg'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  icomm = desc_a%get_mpic()

  call psb_info(ctxt, me, np)

  nglob = desc_a%get_global_rows()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()

  theta = parms%aggr_thresh

  naggr  = nlaggr(me+1)
  ntaggr = sum(nlaggr)

  naggrm1 = sum(nlaggr(1:me))
  naggrp1 = sum(nlaggr(1:me+1))

  filter_mat = (parms%aggr_filter == amg_filter_mat_)

  ! naggr: number of local aggregates
  ! nrow: local rows. 
  ! 
  allocate(adinv(ncol),&
       & omf(ncol),omp(ntaggr),oden(ntaggr),omi(ncol),stat=info)

  if (info /= psb_success_) then 
    info=psb_err_alloc_request_; ierr(1)=6*ncol+ntaggr;
    call psb_errpush(info,name,i_err=ierr,a_err='complex(psb_dpk_)')
    goto 9999      
  end if

  ! Get the diagonal D
  adiag =  a%get_diag(info)
  if (info == psb_success_) &
       & call psb_realloc(ncol,adiag,info)    
  if (info == psb_success_) &
       & call psb_halo(adiag,desc_a,info)
  if (info == psb_success_) call a%cp_to_l(la)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_getdiag')
    goto 9999
  end if

  do i=1,size(adiag)
    if (adiag(i) /= zzero) then
      adinv(i) = zone / adiag(i)
    else
      adinv(i) = zone
    end if
  end do



  ! 1. Allocate Ptilde in sparse matrix form 
  call op_prol%mv_to(tmpcoo)
  call ptilde%mv_from(tmpcoo)
  call ptilde%cscnv(info,type='csr')

  if (info == psb_success_) call la%cscnv(am3,info,type='csr',dupl=psb_dupl_add_)
  if (info == psb_success_) call la%cscnv(da,info,type='csr',dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spcnv')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Initial copies done.'

  call da%scal(adinv,info)

  call psb_spspmm(da,ptilde,dap,info)

  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='spspmm 1')
    goto 9999
  end if

  call dap%clone(atmp,info)

  call psb_sphalo(atmp,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.,outfmt='CSR  ')
  if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=am4)      
  if (info == psb_success_) call am4%free()

  call psb_spspmm(da,atmp,dadap,info)
  call atmp%free()

  !  !$  write(0,*) 'Columns of AP',psb_sp_get_ncols(ap)
  !  !$  write(0,*) 'Columns of ADAP',psb_sp_get_ncols(adap)
  call dap%mv_to(csc_dap)
  call dadap%mv_to(csc_dadap)

  call csc_mat_col_prod(csc_dap,csc_dadap,omp,info)
  call csc_mat_col_prod(csc_dadap,csc_dadap,oden,info)
  call psb_sum(ctxt,omp)
  call psb_sum(ctxt,oden)
  ! !$  write(0,*) trim(name),' OMP :',omp
  ! !$  write(0,*) trim(name),' ODEN:',oden

  omp = omp/oden

  ! !$  write(0,*) 'Check on output prolongator ',omp(1:min(size(omp),10))
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done NUMBMM 1'

  call am3%mv_to(acsr3)
  ! Compute omega_int
  ommx = zzero
  do i=1, ncol
    if (ilaggr(i) >0) then 
      omi(i) = omp(ilaggr(i))
    else
      omi(i) = zzero
    end if
    if(abs(omi(i)) .gt. abs(ommx)) ommx = omi(i)
  end do
  ! Compute omega_fine
  do i=1, nrow
    omf(i) = ommx
    do j=acsr3%irp(i),acsr3%irp(i+1)-1
      if(abs(omi(acsr3%ja(j))) .lt. abs(omf(i))) omf(i)=omi(acsr3%ja(j))
    end do
!!$    if(min(real(omf(i)),aimag(omf(i))) < dzero) omf(i) = zzero
    if(psb_minreal(omf(i)) < dzero) omf(i) = zzero
  end do

  omf(1:nrow) = omf(1:nrow) * adinv(1:nrow)

  if (filter_mat) then
    !
    ! Build the filtered matrix Af from A
    ! 
    call la%cscnv(acsrf,info,dupl=psb_dupl_add_)

    do i=1,nrow
      tmp = zzero
      jd  = -1 
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) jd = j 
        if (abs(acsrf%val(j)) < theta*sqrt(abs(adiag(i)*adiag(acsrf%ja(j))))) then
          tmp=tmp+acsrf%val(j)
          acsrf%val(j)=zzero
        endif
      enddo
      if (jd == -1) then 
        write(0,*) 'Wrong input: we need the diagonal!!!!', i
      else
        acsrf%val(jd)=acsrf%val(jd)-tmp
      end if
    enddo
    ! Take out zeroed terms 
    call acsrf%clean_zeros(info)

    !
    ! Build the smoothed prolongator using the filtered matrix
    ! 
    do i=1,acsrf%get_nrows()
      do j=acsrf%irp(i),acsrf%irp(i+1)-1
        if (acsrf%ja(j) == i) then 
          acsrf%val(j) = zone - omf(i)*acsrf%val(j) 
        else
          acsrf%val(j) = - omf(i)*acsrf%val(j) 
        end if
      end do
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SYMBMM 1'

    call af%mv_from(acsrf)
    !
    ! op_prol = (I-w*D*Af)Ptilde
    ! Doing it this way means to consider diag(Af_i)
    ! 
    !
    call psb_spspmm(af,ptilde,op_prol,info)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done SPSPMM 1'
  else
    !
    ! Build the smoothed prolongator using the original matrix
    !
    do i=1,acsr3%get_nrows()
      do j=acsr3%irp(i),acsr3%irp(i+1)-1
        if (acsr3%ja(j) == i) then 
          acsr3%val(j) = zone - omf(i)*acsr3%val(j) 
        else
          acsr3%val(j) = - omf(i)*acsr3%val(j) 
        end if
      end do
    end do

    call am3%mv_from(acsr3)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done gather, going for SYMBMM 1'
    !
    ! 
    ! op_prol = (I-w*D*A)Ptilde
    ! 
    !
    call psb_spspmm(am3,ptilde,op_prol,info)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Done NUMBMM 1'

  end if


  !
  ! Ok, let's start over with the restrictor
  ! 
  call ptilde%transc(rtilde)
  call la%cscnv(atmp,info,type='csr')
  call psb_sphalo(atmp,desc_a,am4,info,&
       & colcnv=.true.,rowscale=.true.)
  nrt  = am4%get_nrows() 
  call am4%csclip(atmp2,info,lone,nrt,lone,ncol)
  call atmp2%cscnv(info,type='CSR')
  if (info == psb_success_) call psb_rwextd(ncol,atmp,info,b=atmp2)      
  call am4%free()
  call atmp2%free()

  ! This is to compute the transpose. It ONLY works if the
  ! original A has a symmetric pattern.
  call atmp%transc(atmp2) 
  call atmp2%csclip(dat,info,lone,nrow,lone,ncol)
  call dat%cscnv(info,type='csr')
  call dat%scal(adinv,info)

  ! Now for the product. 
  call psb_spspmm(dat,ptilde,datp,info)

  call datp%clone(atmp2,info)
  call psb_sphalo(atmp2,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.,outfmt='CSR  ')
  if (info == psb_success_) call psb_rwextd(ncol,atmp2,info,b=am4)      
  if (info == psb_success_) call am4%free()


  call psb_symbmm(dat,atmp2,datdatp,info)
  call psb_numbmm(dat,atmp2,datdatp)
  call atmp2%free()

  call datp%mv_to(csc_datp)    
  call datdatp%mv_to(csc_datdatp)    

  call csc_mat_col_prod(csc_datp,csc_datdatp,omp,info)
  call csc_mat_col_prod(csc_datdatp,csc_datdatp,oden,info)
  call psb_sum(ctxt,omp)
  call psb_sum(ctxt,oden)


  ! !$  write(debug_unit,*) trim(name),' OMP_R :',omp
  ! ! $  write(debug_unit,*) trim(name),' ODEN_R:',oden
  omp = omp/oden
  ! !$  write(0,*) 'Check on output restrictor',omp(1:min(size(omp),10))
  ! Compute omega_int
  ommx = zzero
  do i=1, ncol
    if (ilaggr(i) >0) then 
      omi(i) = omp(ilaggr(i))
    else
      omi(i) = zzero
    end if
    if(abs(omi(i)) .gt. abs(ommx)) ommx = omi(i)
  end do
  ! Compute omega_fine
  ! Going over the columns of atmp means going over the rows
  ! of A^T. Hopefully ;-) 
  call atmp%cp_to(acsc)

  do i=1, nrow
    omf(i) = ommx
    do j= acsc%icp(i),acsc%icp(i+1)-1
      if(abs(omi(acsc%ia(j))) .lt. abs(omf(i))) omf(i)=omi(acsc%ia(j))
    end do
!!$    if(min(real(omf(i)),aimag(omf(i))) < dzero) omf(i) = zzero
    if(psb_minreal(omf(i)) < dzero) omf(i) = zzero
  end do
  omf(1:nrow) = omf(1:nrow)*adinv(1:nrow)
  call psb_halo(omf,desc_a,info)
  call acsc%free() 


  call atmp%mv_to(acsr1)

  do i=1,acsr1%get_nrows()
    do j=acsr1%irp(i),acsr1%irp(i+1)-1
      if (acsr1%ja(j) == i) then 
        acsr1%val(j) = zone - acsr1%val(j)*omf(acsr1%ja(j))
      else
        acsr1%val(j) =      - acsr1%val(j)*omf(acsr1%ja(j))
      end if
    end do
  end do
  call atmp%mv_from(acsr1)

  call rtilde%mv_to(tmpcoo)
  nzl = tmpcoo%get_nzeros()
  i=0
  do k=1, nzl
    if ((naggrm1 < tmpcoo%ia(k)) .and. (tmpcoo%ia(k) <= naggrp1)) then
      i = i+1
      tmpcoo%val(i) = tmpcoo%val(k)
      tmpcoo%ia(i)  = tmpcoo%ia(k)
      tmpcoo%ja(i)  = tmpcoo%ja(k)
    end if
  end do
  call tmpcoo%set_nzeros(i)
  call rtilde%mv_from(tmpcoo)
  call rtilde%cscnv(info,type='csr')

  call psb_spspmm(rtilde,atmp,op_restr,info)

  !
  ! Now we have to gather the halo of op_prol, and add it to itself
  ! to multiply it by A,
  !
  call op_prol%clone(tmp_prol,info) 
  if (info == psb_success_) call psb_sphalo(tmp_prol,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.)
  if (info == psb_success_) call psb_rwextd(ncol,tmp_prol,info,b=am4)      
  if (info == psb_success_) call am4%free()

  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Halo of op_prol')
    goto 9999
  end if

  !
  ! Now we have to fix this.  The only rows of B that are correct 
  ! are those corresponding to "local" aggregates, i.e. indices in ilaggr(:)
  !
  call op_restr%mv_to(tmpcoo)

  nzl = tmpcoo%get_nzeros()
  i=0
  do k=1, nzl
    if ((naggrm1 < tmpcoo%ia(k)) .and. (tmpcoo%ia(k) <= naggrp1)) then
      i = i+1
      tmpcoo%val(i) = tmpcoo%val(k)
      tmpcoo%ia(i)  = tmpcoo%ia(k)
      tmpcoo%ja(i)  = tmpcoo%ja(k)
    end if
  end do
  call tmpcoo%set_nzeros(i)
  call op_restr%mv_from(tmpcoo)
  call op_restr%cscnv(info,type='csr')


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting sphalo/ rwxtd'

  call psb_spspmm(la,tmp_prol,am3,info)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done SPSPMM 2'

  call psb_sphalo(am3,desc_a,am4,info,&
       & colcnv=.false.,rowscale=.true.)
  if (info == psb_success_) call psb_rwextd(ncol,am3,info,b=am4)      
  if (info == psb_success_) call am4%free()

  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Extend am3')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done sphalo/ rwxtd'

  call psb_spspmm(op_restr,am3,ac,info)
  if (info == psb_success_) call am3%free()
  if (info == psb_success_) call ac%cscnv(info,type='coo',dupl=psb_dupl_add_)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         &a_err='Build ac = op_restr x am3')
    goto 9999
  end if



  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Done smooth_aggregate '
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return


contains

  subroutine csc_mat_col_prod(a,b,v,info)
    implicit none 
    type(psb_lz_csc_sparse_mat), intent(in) :: a, b 
    complex(psb_dpk_), intent(out)             :: v(:)
    integer(psb_ipk_), intent(out)           :: info

    integer(psb_lpk_)                           :: i,j,k, nr, nc,iap,nra,ibp,nrb

    info = psb_success_
    nc   = a%get_ncols()
    if (nc /= b%get_ncols()) then 
      write(0,*) 'Matrices A and B should have same columns'
      info = -1
      return
    end if

    do j=1, nc
      iap  = a%icp(j)
      nra  = a%icp(j+1)-iap
      ibp  = b%icp(j)
      nrb  = b%icp(j+1)-ibp
      v(j) = sparse_srtd_dot(nra,a%ia(iap:iap+nra-1),a%val(iap:iap+nra-1),&
           & nrb,b%ia(ibp:ibp+nrb-1),b%val(ibp:ibp+nrb-1))
    end do

  end subroutine csc_mat_col_prod


  subroutine csr_mat_row_prod(a,b,v,info)
    implicit none 
    type(psb_lz_csr_sparse_mat), intent(in) :: a, b 
    complex(psb_dpk_), intent(out)             :: v(:)
    integer(psb_ipk_), intent(out)           :: info

    integer(psb_lpk_)                        :: i,j,k, nr, nc,iap,nca,ibp,ncb

    info = psb_success_
    nr   = a%get_nrows()
    if (nr /= b%get_nrows()) then 
      write(0,*) 'Matrices A and B should have same rows'
      info = -1
      return
    end if

    do j=1, nr
      iap  = a%irp(j)
      nca  = a%irp(j+1)-iap
      ibp  = b%irp(j)
      ncb  = b%irp(j+1)-ibp
      v(j) = sparse_srtd_dot(nca,a%ja(iap:iap+nca-1),a%val(iap:iap+nca-1),&
           & ncb,b%ja(ibp:ibp+ncb-1),b%val(ibp:ibp+ncb-1))
    end do

  end subroutine csr_mat_row_prod


  function sparse_srtd_dot(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
    implicit none 
    integer(psb_lpk_), intent(in) :: nv1,nv2
    integer(psb_lpk_), intent(in) :: iv1(:), iv2(:)
    complex(psb_dpk_), intent(in) :: v1(:),v2(:)
    complex(psb_dpk_)      :: dot

    integer(psb_lpk_) :: i,j,k, ip1, ip2

    dot = zzero 
    ip1 = 1
    ip2 = 1

    do 
      if (ip1 > nv1) exit
      if (ip2 > nv2) exit
      if (iv1(ip1) == iv2(ip2)) then 
        dot = dot + conjg(v1(ip1))*v2(ip2)
        ip1 = ip1 + 1
        ip2 = ip2 + 1
      else if (iv1(ip1) < iv2(ip2)) then 
        ip1 = ip1 + 1 
      else
        ip2 = ip2 + 1 
      end if
    end do

  end function sparse_srtd_dot

  subroutine local_dump(me,mat,name,header)
    type(psb_lzspmat_type), intent(in) :: mat
    integer(psb_ipk_), intent(in)       :: me
    character(len=*), intent(in)        :: name
    character(len=*), intent(in)        :: header
    character(len=80) :: filename

    write(filename,'(a,a,i0,a,i0,a)') trim(name),'.p',me
    open(20+me,file=filename)
    call mat%print(20+me,head=trim(header))
    close(20+me)
  end subroutine local_dump

end subroutine amg_zaggrmat_minnrg_bld 

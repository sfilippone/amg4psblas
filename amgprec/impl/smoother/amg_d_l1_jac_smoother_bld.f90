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
subroutine amg_d_l1_jac_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)

  use psb_base_mod
  use amg_d_diag_solver
  use amg_d_jac_smoother, amg_protect_name => amg_d_l1_jac_smoother_bld
  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target          :: a
  Type(psb_desc_type), Intent(inout)                 :: desc_a
  class(amg_d_l1_jac_smoother_type), intent(inout)      :: sm
  integer(psb_ipk_), intent(out)                     :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  ! Local variables
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota, nzeros
  real(psb_dpk_), allocatable  :: arwsum(:)
  type(psb_dspmat_type)      :: tmpa
  integer(psb_ipk_) :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20) :: name='d_l1_jac_smoother_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()
  n_col  = desc_a%get_local_cols()
  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()

  if( sm%checkres ) sm%pa => a

  select type (smsv => sm%sv)
  class is (amg_d_diag_solver_type)
    call sm%nd%free()
    sm%pa => a
    sm%nd_nnz_tot = nztota
    
    call psb_sum(ictxt,sm%nd_nnz_tot)
    call sm%sv%build(a,desc_a,info,amold=amold,vmold=vmold)

  class default
    if (smsv%is_global()) then
      ! Do not put anything into SM%ND since the solver
      ! is acting globally.
      call sm%nd%free()
      sm%nd_nnz_tot = 0
      call psb_sum(ictxt,sm%nd_nnz_tot)
      call sm%sv%build(a,desc_a,info,amold=amold,vmold=vmold)
    else

      call a%csclip(tmpa,info,&
           & jmax=nrow_a,rscale=.false.,cscale=.false.)
      
      call a%csclip(sm%nd,info,&
           & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
      
      arwsum = sm%nd%arwsum(info)
      
      call combine_dl1(-done,arwsum,sm%nd,info)
      call combine_dl1(done,arwsum,tmpa,info)
      
      sm%nd_nnz_tot = sm%nd%get_nzeros()
      call psb_sum(ictxt,sm%nd_nnz_tot)

      call sm%sv%build(tmpa,desc_a,info,amold=amold,vmold=vmold)

      if (info == psb_success_) then
        if (present(amold)) then
          call sm%nd%cscnv(info,&
               & mold=amold,dupl=psb_dupl_add_)
        else
          call sm%nd%cscnv(info,&
               & type='csr',dupl=psb_dupl_add_)
        endif
      end if
    end if
  end select
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if


  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='solver build')
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
contains
  
  subroutine combine_dl1(alpha,dl1,mat,info)
    implicit none
    real(psb_dpk_), intent(in) ::  alpha, dl1(:)
    type(psb_dspmat_type), intent(inout) :: mat
    integer(psb_ipk_), intent(out) :: info 
    !
    integer(psb_ipk_) :: k, nz, nrm, dp
    type(psb_d_coo_sparse_mat) :: tcoo

    call mat%mv_to(tcoo)
    nz = tcoo%get_nzeros()
    nrm = min(size(dl1),tcoo%get_nrows(),tcoo%get_ncols())
!!$    write(0,*) 'Check on combine_dl1: ',nrm, tcoo%get_nrows(),tcoo%get_ncols(), nz
    call tcoo%ensure_size(nz+nrm)
    call tcoo%set_dupl(psb_dupl_add_)
    do k=1,nrm
      if (dl1(k) /= dzero) then
        nz = nz + 1
        tcoo%ia(nz)  = k
        tcoo%ja(nz)  = k
        tcoo%val(nz) = alpha*dl1(k)
      end if
    end do
    call tcoo%set_nzeros(nz)
    call tcoo%fix(info)
    call mat%mv_from(tcoo)
  end subroutine combine_dl1


end subroutine amg_d_l1_jac_smoother_bld

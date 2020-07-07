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
subroutine mld_s_as_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)
  
  use psb_base_mod
  use mld_s_as_smoother, mld_protect_nam => mld_s_as_smoother_bld
  Implicit None

  ! Arguments
  type(psb_sspmat_type), intent(in), target        :: a
  Type(psb_desc_type), Intent(inout)                 :: desc_a 
  class(mld_s_as_smoother_type), intent(inout)     :: sm
  integer(psb_ipk_), intent(out)                     :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  ! Local variables
  type(psb_sspmat_type) :: blck, atmp
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nhalo, novr, data_
  real(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer(psb_ipk_) :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20) :: name='s_as_smoother_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  sm%pa  => a
  novr   = sm%novr
  if (novr < 0) then
    info=psb_err_invalid_ovr_num_
    call psb_errpush(info,name,&
         & i_err=(/novr,izero,izero,izero,izero,izero/))
    goto 9999
  endif

  if ((novr == 0).or.(np == 1)) then 
    call psb_cdcpy(desc_a,sm%desc_data,info)
    If(debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & '  done cdcpy'
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_cdcpy'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Early return: P>=3 N_OVR=0'
    call blck%csall(izero,izero,info,ione)
  else

    !
    ! Build the auxiliary descriptor desc_p%matrix_data(psb_n_row_).
    ! This is done by psb_cdbldext (interface to psb_cdovr), which is
    ! independent of CSR, and has been placed in the tools directory
    ! of PSBLAS, instead of the mlprec directory of MLD2P4, because it
    ! might be used independently of the AS preconditioner, to build
    ! a descriptor for an extended stencil in a PDE solver. 
    !
    call psb_cdbldext(a,desc_a,novr,sm%desc_data,info,extype=psb_ovt_asov_)
    if(debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' From cdbldext _:',sm%desc_data%get_local_rows(),&
         & sm%desc_data%get_local_cols()
    
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_cdbldext'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'Before sphalo '
    
    !
    ! Retrieve the remote sparse matrix rows required for the AS extended
    ! matrix
    data_ = psb_comm_ext_
    Call psb_sphalo(a,sm%desc_data,blck,info,data=data_,rowscale=.true.)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sphalo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    
    if (debug_level >=psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'After psb_sphalo ',&
         & blck%get_nrows(), blck%get_nzeros()

  End if
  if (info == psb_success_) &
       & call sm%sv%build(a,sm%desc_data,info,&
       &  blck,amold=amold,vmold=vmold)

  nrow_a = a%get_nrows()
  n_row  = sm%desc_data%get_local_rows()
  n_col  = sm%desc_data%get_local_cols()

  if (info == psb_success_) call a%csclip(sm%nd,info,&
       & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
  if (info == psb_success_) call blck%csclip(atmp,info,&
       & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
  if (info == psb_success_) call psb_rwextd(n_row,sm%nd,info,b=atmp) 

  if (info == psb_success_) then 
    if (present(amold)) then 
      call sm%nd%cscnv(info,&
           & mold=amold,dupl=psb_dupl_add_)
    else
      call sm%nd%cscnv(info,&
           & type='csr',dupl=psb_dupl_add_)
    end if
  end if
  if (info == psb_success_) then 
    if (present(imold)) then 
      call sm%desc_data%cnv(imold)
    end if
  end if
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if
  sm%nd_nnz_tot = sm%nd%get_nzeros()
  call psb_sum(ictxt,sm%nd_nnz_tot)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine mld_s_as_smoother_bld

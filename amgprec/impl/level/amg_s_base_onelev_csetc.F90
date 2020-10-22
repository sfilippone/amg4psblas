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
subroutine amg_s_base_onelev_csetc(lv,what,val,info,pos,idx)
  
  use psb_base_mod
  use amg_s_onelev_mod, amg_protect_name => amg_s_base_onelev_csetc
  use amg_s_base_aggregator_mod
  use amg_s_dec_aggregator_mod
  use amg_s_symdec_aggregator_mod
  use amg_s_jac_smoother
  use amg_s_as_smoother
  use amg_s_diag_solver
  use amg_s_l1_diag_solver
  use amg_s_ilu_solver
  use amg_s_id_solver
  use amg_s_gs_solver
#if defined(HAVE_SLU_)
  use amg_s_slu_solver
#endif
#if defined(HAVE_MUMPS_)
  use amg_s_mumps_solver
#endif

  Implicit None

  ! Arguments
  class(amg_s_onelev_type), intent(inout) :: lv 
  character(len=*), intent(in)              :: what 
  character(len=*), intent(in)              :: val
  integer(psb_ipk_), intent(out)            :: info
  character(len=*), optional, intent(in)    :: pos
  integer(psb_ipk_), intent(in), optional   :: idx
  ! Local 
  integer(psb_ipk_)  :: ipos_, err_act
  character(len=20) :: name='s_base_onelev_csetc'
  integer(psb_ipk_) :: ival 
  type(amg_s_base_smoother_type)   :: amg_s_base_smoother_mold
  type(amg_s_jac_smoother_type)    ::  amg_s_jac_smoother_mold
  type(amg_s_l1_jac_smoother_type) ::  amg_s_l1_jac_smoother_mold
  type(amg_s_as_smoother_type)     ::  amg_s_as_smoother_mold
  type(amg_s_diag_solver_type)     ::  amg_s_diag_solver_mold
  type(amg_s_l1_diag_solver_type)  ::  amg_s_l1_diag_solver_mold
  type(amg_s_ilu_solver_type)      ::  amg_s_ilu_solver_mold
  type(amg_s_id_solver_type)       ::  amg_s_id_solver_mold
  type(amg_s_gs_solver_type)       ::  amg_s_gs_solver_mold
  type(amg_s_bwgs_solver_type)     ::  amg_s_bwgs_solver_mold
#if defined(HAVE_SLU_)
  type(amg_s_slu_solver_type)   ::  amg_s_slu_solver_mold
#endif
#if defined(HAVE_MUMPS_)
  type(amg_s_mumps_solver_type) ::  amg_s_mumps_solver_mold
#endif
  

  call psb_erractionsave(err_act)

  info = psb_success_

  ival = lv%stringval(val)


  if (present(pos)) then
    select case(psb_toupper(trim(pos)))
    case('PRE')
      ipos_ = amg_smooth_pre_
    case('POST')
      ipos_ = amg_smooth_post_
    case default
      ipos_ = amg_smooth_both_
    end select
  else
    ipos_ = amg_smooth_both_
  end if
  
  select case (psb_toupper(trim(what)))
  case ('SMOOTHER_TYPE')
    select case (psb_toupper(trim(val)))
    case ('NOPREC','NONE')
      call lv%set(amg_s_base_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(amg_s_id_solver_mold,info,pos=pos)
      
    case ('JAC','JACOBI')
      call lv%set(amg_s_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(amg_s_diag_solver_mold,info,pos=pos)

    case ('L1-JACOBI')
      call lv%set(amg_s_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(amg_s_l1_diag_solver_mold,info,pos=pos)
      
    case ('BJAC')
      call lv%set(amg_s_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(amg_s_ilu_solver_mold,info,pos=pos)
      
    case ('L1-BJAC')
      call lv%set(amg_s_l1_jac_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(amg_s_ilu_solver_mold,info,pos=pos)

    case ('AS')
      call lv%set(amg_s_as_smoother_mold,info,pos=pos)
      if (info == 0) call lv%set(amg_s_ilu_solver_mold,info,pos=pos)

    case ('GS','FWGS')
      call lv%set(amg_s_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(amg_s_gs_solver_mold,info,pos='pre')
      if (allocated(lv%sm2a)) deallocate(lv%sm2a)
    case ('BWGS')
      call lv%set(amg_s_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(amg_s_bwgs_solver_mold,info,pos='pre')
      if (allocated(lv%sm2a)) deallocate(lv%sm2a)
    case ('FBGS')
      call lv%set(amg_s_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(amg_s_gs_solver_mold,info,pos='pre')
      call lv%set(amg_s_jac_smoother_mold,info,pos='post')
      if (info == 0) call lv%set(amg_s_bwgs_solver_mold,info,pos='post')
    case ('L1-GS','L1-FWGS')
      call lv%set(amg_s_l1_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(amg_s_gs_solver_mold,info,pos='pre')
      if (allocated(lv%sm2a)) deallocate(lv%sm2a)
    case ('L1-BWGS')
      call lv%set(amg_s_l1_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(amg_s_bwgs_solver_mold,info,pos='pre')
      if (allocated(lv%sm2a)) deallocate(lv%sm2a)      
    case ('L1-FBGS')
      call lv%set(amg_s_l1_jac_smoother_mold,info,pos='pre')
      if (info == 0) call lv%set(amg_s_gs_solver_mold,info,pos='pre')
      call lv%set(amg_s_l1_jac_smoother_mold,info,pos='post')
      if (info == 0) call lv%set(amg_s_bwgs_solver_mold,info,pos='post')
      
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    if ((ipos_==amg_smooth_pre_).or.(ipos_==amg_smooth_both_)) then 
      if (allocated(lv%sm)) call lv%sm%default()
    end if
    if ((ipos_==amg_smooth_post_).or.(ipos_==amg_smooth_both_)) then
      if (allocated(lv%sm2a)) call lv%sm2a%default()
    end if
    

  case('SUB_SOLVE')
    select case (psb_toupper(trim(val)))
    case ('NONE','NOPREC','FACT_NONE')
      call lv%set(amg_s_id_solver_mold,info,pos=pos)
      
    case ('DIAG')
      call lv%set(amg_s_diag_solver_mold,info,pos=pos)
      
    case ('L1-DIAG')
      call lv%set(amg_s_l1_diag_solver_mold,info,pos=pos)
      
    case ('GS','FGS','FWGS')
      call lv%set(amg_s_gs_solver_mold,info,pos=pos)
      
    case ('BGS','BWGS')
      call lv%set(amg_s_bwgs_solver_mold,info,pos=pos)
      
    case ('ILU','ILUT','MILU')
      call lv%set(amg_s_ilu_solver_mold,info,pos=pos)
      if (info == 0) then
        if ((ipos_==amg_smooth_pre_) .or.(ipos_==amg_smooth_both_)) then 
          call lv%sm%sv%set('SUB_SOLVE',val,info)
        end if
        if ((ipos_==amg_smooth_post_).or.(ipos_==amg_smooth_both_))then 
          if (allocated(lv%sm2a)) call lv%sm2a%sv%set('SUB_SOLVE',val,info)
        end if
      end if
#ifdef HAVE_SLU_
    case ('SLU') 
      call lv%set(amg_s_slu_solver_mold,info,pos=pos)
#endif
#ifdef HAVE_MUMPS_
    case ('MUMPS') 
      call lv%set(amg_s_mumps_solver_mold,info,pos=pos)
#endif
    case default
      !
      ! Do nothing and hope for the best :) 
      !
    end select
    
  case ('ML_CYCLE')
    lv%parms%ml_cycle      = amg_stringval(val)

  case ('PAR_AGGR_ALG')
    ival = amg_stringval(val)
    lv%parms%par_aggr_alg  = ival
    if (allocated(lv%aggr)) then
      call lv%aggr%free(info)
      if (info == 0) deallocate(lv%aggr,stat=info)
      if (info /= 0) then
        info = psb_err_internal_error_
        return
      end if
    end if
    
    select case(ival)
    case(amg_dec_aggr_)
      allocate(amg_s_dec_aggregator_type :: lv%aggr, stat=info)
    case(amg_sym_dec_aggr_)
      allocate(amg_s_symdec_aggregator_type :: lv%aggr, stat=info)
    case default
      info =  psb_err_internal_error_
    end select
    if (info == psb_success_) call lv%aggr%default()
    
  case ('AGGR_ORD')
    lv%parms%aggr_ord      = amg_stringval(val)

  case ('AGGR_TYPE')
    lv%parms%aggr_type     = amg_stringval(val)
    if (allocated(lv%aggr)) call lv%aggr%set_aggr_type(lv%parms,info)

  case ('AGGR_PROL')
    lv%parms%aggr_prol     = amg_stringval(val)

  case ('COARSE_MAT')
    lv%parms%coarse_mat    = amg_stringval(val)

  case ('AGGR_OMEGA_ALG')
    lv%parms%aggr_omega_alg= amg_stringval(val)

  case ('AGGR_EIG')
    lv%parms%aggr_eig      = amg_stringval(val)

  case ('AGGR_FILTER')
    lv%parms%aggr_filter   = amg_stringval(val)

  case ('COARSE_SOLVE')
    lv%parms%coarse_solve    = amg_stringval(val)

  case default
    if ((ipos_==amg_smooth_pre_) .or.(ipos_==amg_smooth_both_)) then 
      if (allocated(lv%sm)) then 
        call lv%sm%set(what,val,info,idx=idx)
      end if
    end if
    if ((ipos_==amg_smooth_post_).or.(ipos_==amg_smooth_both_))then 
      if (allocated(lv%sm2a)) then 
        call lv%sm2a%set(what,val,info,idx=idx)
      end if
    end if
    if (allocated(lv%aggr)) call lv%aggr%set(what,val,info,idx=idx)

  end select


  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_s_base_onelev_csetc

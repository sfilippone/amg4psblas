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
! File: amg_dprecset.f90
!
! Subroutine: amg_dprecseti
! Version: real
!
!  This routine sets the integer parameters defining the preconditioner. More
!  precisely, the integer parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level.
!
!  To set character and real parameters, see amg_dprecsetc and amg_dprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(amg_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the AMG4PSBLAS User's and Reference Guide.
!    val     -  integer, input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in the AMG4PSBLAS User's and Reference Guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set.
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  AMG4PSBLAS developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface amg_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see amg_prec_mod.f90).
!
subroutine amg_dcprecseti(p,what,val,info,ilev,ilmax,pos,idx)

  use psb_base_mod
  use amg_d_prec_mod, amg_protect_name => amg_dcprecseti
  use amg_d_jac_smoother
  use amg_d_as_smoother
  use amg_d_diag_solver
  use amg_d_l1_diag_solver
  use amg_d_ilu_solver
  use amg_d_id_solver
  use amg_d_gs_solver
  use amg_d_ainv_solver
  use amg_d_invk_solver
  use amg_d_invt_solver
#if defined(HAVE_UMF_)
  use amg_d_umf_solver
#endif
#if defined(HAVE_SLUDIST_)
  use amg_d_sludist_solver
#endif
#if defined(HAVE_SLU_)
  use amg_d_slu_solver
#endif
#if defined(HAVE_MUMPS_)
  use amg_d_mumps_solver
#endif


  implicit none

  ! Arguments
  class(amg_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)            :: what
  integer(psb_ipk_), intent(in)           :: val
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), optional, intent(in) :: ilev,ilmax
  character(len=*), optional, intent(in)  :: pos
  integer(psb_ipk_), intent(in), optional :: idx

  ! Local variables
  integer(psb_ipk_)                      :: ilev_, nlev_, ilmax_, il
  character(len=*), parameter            :: name='amg_precseti'

  info = psb_success_

  if (.not.allocated(p%precv)) then
    info = 3111
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call amg_PRECINIT'
    return
  endif

  nlev_ = size(p%precv)

  if (present(ilev)) then
    ilev_ = ilev
    if (present(ilmax)) then
      ilmax_ = ilmax
    else
      ilmax_ = ilev_
    end if
  else
    ilev_  = 1
    ilmax_ = nlev_
  end if
  if ((ilev_<1).or.(ilev_ > nlev_)) then
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if ((ilmax_<1).or.(ilmax_ > nlev_)) then
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILMAX/NLEV combination',ilmax_, nlev_
    return
  endif

  select case(psb_toupper(what))
  case ('MIN_COARSE_SIZE')
    p%ag_data%min_coarse_size = max(val,-1)
    return
  case ('MIN_COARSE_SIZE_PER_PROCESS')
    p%ag_data%min_coarse_size_per_process = max(val,-1)
    return
  case('MAX_LEVS')
    p%ag_data%max_levs = max(val,1)
    return
  case ('OUTER_SWEEPS')
    p%outer_sweeps = max(val,1)
    return
  end select

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then

      select case(psb_toupper(what))
      case('SUB_OVR','SUB_FILLIN','SMOOTHER_SWEEPS')
        do il=ilev_, ilmax_
          call p%precv(il)%set(what,val,info,pos=pos)
        end do


      case('COARSE_SWEEPS')
        if (ilev_ /= nlev_) then
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(nlev_)%set('SMOOTHER_SWEEPS',val,info,pos=pos)

      case('COARSE_FILLIN')
        if (ilev_ /= nlev_) then
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(nlev_)%set('SUB_FILLIN',val,info,pos=pos)

     case('BJAC_ITRACE')
       if (ilev_ /= nlev_) then
         write(psb_err_unit,*) name,&
              & ': Error: Inconsistent specification of WHAT vs. ILEV'
         info = -2
         return
       end if
       call p%precv(nlev_)%set('SMOOTHER_ITRACE',val,info,pos=pos)

    case('BJAC_RESCHECK')
      if (ilev_ /= nlev_) then
        write(psb_err_unit,*) name,&
             & ': Error: Inconsistent specification of WHAT vs. ILEV'
        info = -2
        return
      end if
      call p%precv(nlev_)%set('SMOOTHER_RESIDUAL',val,info,pos=pos)

      case default
        do il=ilev_, ilmax_
          call p%precv(il)%set(what,val,info,pos=pos,idx=idx)
        end do
      end select


  else if (.not.present(ilev)) then
    !
    ! ilev not specified: set preconditioner parameters at all the appropriate
    ! levels
    !
    select case(psb_toupper(trim(what)))
    case('SUB_OVR','SUB_FILLIN','SMOOTHER_SWEEPS')
      do ilev_=1,max(1,nlev_-1)
        call p%precv(ilev_)%set(what,val,info,pos=pos)
        if (info /= 0) return
      end do


    case('COARSE_SWEEPS')

      if (nlev_ > 1) then
        call p%precv(nlev_)%set('SMOOTHER_SWEEPS',val,info,pos=pos)
      end if

    case('COARSE_FILLIN')
      if (nlev_ > 1) then
        call p%precv(nlev_)%set('SUB_FILLIN',val,info,pos=pos)
      end if

    case('BJAC_ITRACE')
      if (nlev_ > 1) then
        call p%precv(nlev_)%set('SMOOTHER_ITRACE',val,info,pos=pos)
      end if
    case('BJAC_RESCHECK')
      if (nlev_ > 1) then
        call p%precv(nlev_)%set('SMOOTHER_RESIDUAL',val,info,pos=pos)
      end if
    case default
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,val,info,pos=pos,idx=idx)
      end do
    end select

  endif

end subroutine amg_dcprecseti

!
! Subroutine: amg_dprecsetc
! Version: real
!
!  This routine sets the character parameters defining the preconditioner. More
!  precisely, the character parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level.
!
!  To set integer and real parameters, see amg_dprecseti and amg_dprecsetr,
!  respectively.
!
!
! Arguments:
!    p       -  type(amg_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the AMG4PSBLAS User's and Reference Guide.
!    string  -  character(len=*), input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in the AMG4PSBLAS User's and Reference Guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set.
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  AMG4PSBLAS developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface amg_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see amg_prec_mod.f90).
!
subroutine amg_dcprecsetc(p,what,string,info,ilev,ilmax,pos,idx)

  use psb_base_mod
  use amg_d_prec_mod, amg_protect_name => amg_dcprecsetc
  use amg_d_jac_smoother
  use amg_d_as_smoother
  use amg_d_diag_solver
  use amg_d_l1_diag_solver
  use amg_d_ilu_solver
  use amg_d_id_solver
  use amg_d_gs_solver
  use amg_d_ainv_solver
  use amg_d_invk_solver
  use amg_d_invt_solver
#if defined(HAVE_UMF_)
  use amg_d_umf_solver
#endif
#if defined(HAVE_SLUDIST_)
  use amg_d_sludist_solver
#endif
#if defined(HAVE_SLU_)
  use amg_d_slu_solver
#endif
#if defined(HAVE_MUMPS_)
  use amg_d_mumps_solver
#endif
  use amg_d_krm_solver, only : amg_d_krm_solver_type


  implicit none

  ! Arguments
  class(amg_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)            :: what
  character(len=*), intent(in)            :: string
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), optional, intent(in) :: ilev,ilmax
  character(len=*), optional, intent(in)   :: pos
  integer(psb_ipk_), intent(in), optional  :: idx

  ! Local variables
  integer(psb_ipk_)                      :: ilev_, nlev_,val,ilmax_, il
  character(len=*), parameter            :: name='amg_precsetc'

  info = psb_success_

  if (.not.allocated(p%precv)) then
    info = 3111
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call amg_PRECINIT'
    return
  endif

  nlev_ = size(p%precv)

  if (present(ilev)) then
    ilev_ = ilev
    if (present(ilmax)) then
      ilmax_ = ilmax
    else
      ilmax_ = ilev_
    end if
  else
    ilev_  = 1
    ilmax_ = nlev_
  end if
  if ((ilev_<1).or.(ilev_ > nlev_)) then
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if ((ilmax_<1).or.(ilmax_ > nlev_)) then
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILMAX/NLEV combination',ilmax_, nlev_
    return
  endif

  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then

      select case(psb_toupper(what))
      case('SMOOTHER_TYPE')
        ! Select the type of smoother between the one implemented in the library
        ! every new smoother should be added here
        select case(psb_toupper(string))
        case ('BJAC')
          do il=ilev_, ilmax_
            call  p%precv(il)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          end do
        case ('L1-BJAC')
          do il=ilev_, ilmax_
            call  p%precv(il)%set('SMOOTHER_TYPE',amg_l1_bjac_,info,pos=pos)
          end do
        case('GS','FWGS','FBGS')
          do il=ilev_, ilmax_
            call p%precv(il)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(il)%set('SUB_SOLVE',amg_gs_,info,pos=pos)
          end do
        case('L1-GS','L1-FWGS','L1-FBGS')
          do il=ilev_, ilmax_
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_l1_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_gs_,info,pos=pos)
          end do
        case ('BWGS')
          do il=ilev_, ilmax_
            call p%precv(il)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(il)%set('SUB_SOLVE',amg_bwgs_,info,pos=pos)
          end do
        case('JACOBI')
          do il=ilev_, ilmax_
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_diag_scale_,info,pos=pos)
          end do
        case('L1-JACOBI')
          do il=ilev_, ilmax_
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_l1_diag_scale_,info,pos=pos)
          end do
        end select
      case('SUB_SOLVE','ML_CYCLE','PAR_AGGR_ALG','AGGR_TYPE','SUB_RESTR'&
        & ,'SUB_PROL')
        ! These are handled elsewhere
        do il=ilev_, ilmax_
          call p%precv(nlev_)%set(what,string,info,pos=pos)
        end do
      case('COARSE_MAT')
        ! Select if the coarsest matrix is handled in a distributed way, few
        ! rows per rank, or if it is replicated completely on every rank
        select case(psb_toupper(string))
        case('DISTR')
          do il=ilev_, ilmax_
            call p%precv(il)%set(what,amg_distr_mat_,info,pos=pos)
          end do
        case('REPL')
          do il=ilev_, ilmax_
            call p%precv(il)%set(what,amg_repl_mat_,info,pos=pos)
          end do
      case('BJAC_STOP')
        do il=ilev_, ilmax_
          call p%precv(il)%set('SMOOTHER_STOP',string,info,pos=pos)
        end do
      case('BJAC_TRACE')
        do il=ilev_, ilmax_
          call p%precv(il)%set('SMOOTHER_TRACE',string,info,pos=pos)
        end do
      case('COARSE_SUBSOLVE')
        if (ilev_ /= nlev_) then
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        call p%precv(ilev_)%set('SUB_SOLVE',string,info,pos=pos)
      case('COARSE_SOLVE')
        if (ilev_ /= nlev_) then
          write(psb_err_unit,*) name,&
               & ': Error: Inconsistent specification of WHAT vs. ILEV'
          info = -2
          return
        end if
        if (nlev_ > 1) then
          select case (psb_toupper(string))
          case('BJAC')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          case('L1-BJAC')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_l1_bjac_,info,pos=pos)
#if defined(HAVE_UMF_)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_umf_,info,pos=pos)
#elif defined(HAVE_SLU_)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
#else
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
#endif
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info)
          case('SLU')
#if defined(HAVE_SLU_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif
          case('ILU')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
          case('ILUT')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_t_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
          case('MILU')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_milu_n_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
          case('MUMPS')
#if defined(HAVE_MUMPS_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif
          case('UMF')
#if defined(HAVE_UMF_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_umf_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_SLU_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif

          case('SLUDIST')
#if defined(HAVE_SLUDIST_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_sludist_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#elif defined(HAVE_UMF_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_umf_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_SLU_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif
          case('JACOBI','JAC')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_diag_scale_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)

          case('L1-JACOBI','L1-JAC')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_l1_diag_scale_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
          case('GS','FBGS')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_gs_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
          case('BWGS')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_bwgs_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
          case('L1-GS','L1-FBGS')
            call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_l1_bjac_,info,pos=pos)
            call p%precv(nlev_)%set('SUB_SOLVE',amg_gs_,info,pos=pos)
            call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
          case('KRM')
            block
              type(amg_d_krm_solver_type)          :: krm_slv
              call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_krm_,info,pos=pos)
              call p%precv(nlev_)%set(krm_slv,info)
              call p%precv(nlev_)%default()
              call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
            end block
          end select

      end if

      end select

      case default
        do il=ilev_, ilmax_
          call p%precv(il)%set(what,string,info,pos=pos,idx=idx)
        end do
      end select




  else if (.not.present(ilev)) then
    !
    ! ilev not specified: set preconditioner parameters at all the appropriate
    ! levels
    !
    select case(psb_toupper(trim(what)))
    case('SUB_SOLVE','SUB_RESTR','SUB_PROL',&
         & 'SMOOTHER_TYPE')
      do ilev_=1,max(1,nlev_-1)
        call p%precv(ilev_)%set(what,string,info,pos=pos)
        if (info /= 0) return
      end do

    case('ML_CYCLE','PAR_AGGR_ALG','AGGR_ORD','AGGR_PROL','AGGR_TYPE',&
         & 'AGGR_OMEGA_ALG','AGGR_EIG','AGGR_FILTER')
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,string,info,pos=pos)
        if (info /= 0) return
      end do

    case('COARSE_MAT')
      if (nlev_ > 1) then
        call p%precv(nlev_)%set('COARSE_MAT',string,info,pos=pos)
      end if

    case('COARSE_SOLVE')
      if (nlev_ > 1) then
        call p%precv(nlev_)%set('COARSE_SOLVE',string,info,pos=pos)
        select case (psb_toupper(trim(string)))
        case('BJAC')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
        case('L1-BJAC')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_l1_bjac_,info,pos=pos)
#if defined(HAVE_UMF_)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_umf_,info,pos=pos)
#elif defined(HAVE_SLU_)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
#else
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
#endif
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info)
        case('SLU')
#if defined(HAVE_SLU_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif
        case('ILU')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
        case('ILUT')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_t_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
        case('MILU')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_milu_n_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
        case('MUMPS')
#if defined(HAVE_MUMPS_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif
        case('UMF')
#if defined(HAVE_UMF_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_umf_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_SLU_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif

        case('SLUDIST')
#if defined(HAVE_SLUDIST_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_sludist_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#elif defined(HAVE_UMF_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_umf_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_SLU_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_slu_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_repl_mat_,info,pos=pos)
#elif defined(HAVE_MUMPS_)
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_mumps_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#else
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',psb_ilu_n_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
#endif
        case('JACOBI','JAC')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_diag_scale_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)

        case('L1-JACOBI','L1-JAC')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_l1_diag_scale_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
        case('GS','FBGS')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_gs_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
        case('BWGS')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_bwgs_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
        case('L1-FBGS','L1-GS')
          call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_l1_bjac_,info,pos=pos)
          call p%precv(nlev_)%set('SUB_SOLVE',amg_gs_,info,pos=pos)
          call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
          case('KRM')
            block
              type(amg_d_krm_solver_type)          :: krm_slv
              call p%precv(nlev_)%set('SMOOTHER_TYPE',amg_krm_,info,pos=pos)
              call p%precv(nlev_)%set(krm_slv,info)
              call p%precv(nlev_)%default()
              call p%precv(nlev_)%set('COARSE_MAT',amg_distr_mat_,info,pos=pos)
            end block
        end select
    endif

    case('BJAC_STOP')
      if (nlev_ > 1) then
        call p%precv(il)%set('SMOOTHER_STOP',string,info,pos=pos)
      end if
    case('BJAC_TRACE')
      if (nlev_ > 1) then
        call p%precv(il)%set('SMOOTHER_TRACE',string,info,pos=pos)
      end if

    case('COARSE_SUBSOLVE')
      if (nlev_ > 1) then
        call p%precv(nlev_)%set('SUB_SOLVE',string,info,pos=pos)
      endif

    case default
      do ilev_=1,nlev_
        call p%precv(ilev_)%set(what,string,info,pos=pos,idx=idx)
      end do
    end select

  endif


end subroutine amg_dcprecsetc


!
! Subroutine: amg_dprecsetr
! Version: real
!
!  This routine sets the real parameters defining the preconditioner. More
!  precisely, the real parameter identified by 'what' is assigned the value
!  contained in 'val'.
!  For the multilevel preconditioners, the levels are numbered in increasing
!  order starting from the finest one, i.e. level 1 is the finest level.
!
!  To set integer and character parameters, see amg_dprecseti and amg_dprecsetc,
!  respectively.
!
! Arguments:
!    p       -  type(amg_dprec_type), input/output.
!               The preconditioner data structure.
!    what    -  integer, input.
!               The number identifying the parameter to be set.
!               A mnemonic constant has been associated to each of these
!               numbers, as reported in the AMG4PSBLAS User's and Reference Guide.
!    val     -  real(psb_dpk_), input.
!               The value of the parameter to be set. The list of allowed
!               values is reported in the AMG4PSBLAS User's and Reference Guide.
!    info    -  integer, output.
!               Error code.
!    ilev    -  integer, optional, input.
!               For the multilevel preconditioner, the level at which the
!               preconditioner parameter has to be set.
!               If nlev is not present, the parameter identified by 'what'
!               is set at all the appropriate levels.
!
!  NOTE: currently, the use of the argument ilev is not "safe" and is reserved to
!  AMG4PSBLAS developers. Indeed, by using ilev it is possible to set different values
!  of the same parameter at different levels 1,...,nlev-1, even in cases where
!  the parameter must have the same value at all the levels but the coarsest one.
!  For this reason, the interface amg_precset to this routine has been built in
!  such a way that ilev is not visible to the user (see amg_prec_mod.f90).
!
subroutine amg_dcprecsetr(p,what,val,info,ilev,ilmax,pos,idx)

  use psb_base_mod
  use amg_d_prec_mod, amg_protect_name => amg_dcprecsetr

  implicit none

  ! Arguments
  class(amg_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)            :: what
  real(psb_dpk_), intent(in)              :: val
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), optional, intent(in) :: ilev,ilmax
  character(len=*), optional, intent(in)   :: pos
  integer(psb_ipk_), intent(in), optional  :: idx

  ! Local variables
  integer(psb_ipk_)                      :: ilev_,nlev_, ilmax_, il
  real(psb_dpk_)                         :: thr
  character(len=*), parameter            :: name='amg_precsetr'

  info = psb_success_

  if (present(ilev)) then
    ilev_ = ilev
  else
    ilev_ = 1
  end if

  select case(psb_toupper(what))
  case ('MIN_CR_RATIO')
    p%ag_data%min_cr_ratio = max(done,val)
    return
  end select

  if (.not.allocated(p%precv)) then
    write(psb_err_unit,*) name,&
         &': Error: uninitialized preconditioner,',&
         &' should call amg_PRECINIT'
    info = 3111
    return
  endif
  nlev_ = size(p%precv)

  if (present(ilev)) then
    ilev_ = ilev
    if (present(ilmax)) then
      ilmax_ = ilmax
    else
      ilmax_ = ilev_
    end if
  else
    ilev_  = 1
    ilmax_ = nlev_
  end if
  if ((ilev_<1).or.(ilev_ > nlev_)) then
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if ((ilmax_<1).or.(ilmax_ > nlev_)) then
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILMAX/NLEV combination',ilmax_, nlev_
    return
  endif


  !
  ! Set preconditioner parameters at level ilev.
  !
  if (present(ilev)) then

    select case(psb_toupper(trim(what)))
    case('BJAC_STOPTOL')
      do il=ilev_, ilmax_
        call p%precv(il)%set('SMOOTHER_STOPTOL',val,info,pos=pos)
      end do
    case default
      ! Nothing to do here, setted elsewhere
      do il=ilev_, ilmax_
        call p%precv(il)%set(what,val,info,pos=pos,idx=idx)
      end do
    end select
  else if (.not.present(ilev)) then
    !
    ! ilev not specified: set preconditioner parameters at all the appropriate levels
    !

    select case(psb_toupper(what))
    case('COARSE_ILUTHRS')
      ilev_=nlev_
      call p%precv(ilev_)%set('SUB_ILUTHRS',val,info,pos=pos)
    case('BJAC_STOPTOL')
      ilev_=nlev_
      call p%precv(ilev_)%set('SMOOTHER_STOPTOL',val,info,pos=pos)
    case default

      do il=1,nlev_
        call p%precv(il)%set(what,val,info,pos=pos,idx=idx)
      end do
    end select

  endif

end subroutine amg_dcprecsetr

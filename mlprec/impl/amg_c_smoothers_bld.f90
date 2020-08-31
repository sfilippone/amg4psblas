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
! File: amg_c_smoothers_bld.f90
!
! Subroutine: amg_c_smoothers_bld
! Version:    complex
!
!  This routine performs the final phase of the multilevel preconditioner
!  build process: builds the "smoother" objects at each level,
!  based on the matrix hierarchy prepared by amg_c_hierarchy_bld.
!  
!  A multilevel preconditioner is regarded as an array of 'one-level'
!  data structures, each containing the part of the
!  preconditioner associated to a certain level,
!  (for more details see the description of amg_Tonelev_type in amg_prec_type.f90).
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level. No transfer operators are associated to level 1.
!  Each level provides a "build" method; for the base type, the "one-level"
!  build procedure simply invokes the build method of the first smoother object,
!  and also on the second object if allocated. 
! 
!
! Arguments:
!    a       -  type(psb_cspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_cprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!
!    amold   -  class(psb_c_base_sparse_mat), input, optional
!               Mold for the inner format of matrices contained in the
!               preconditioner
!
!
!    vmold   -  class(psb_c_base_vect_type), input, optional
!               Mold for the inner format of vectors contained in the
!               preconditioner
!
!
!  
subroutine amg_c_smoothers_bld(a,desc_a,prec,info,amold,vmold,imold)

  use psb_base_mod
  !use amg_c_inner_mod
  use amg_c_prec_mod, amg_protect_name => amg_c_smoothers_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  class(amg_cprec_type),intent(inout),target         :: prec
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_c_base_sparse_mat), intent(in), optional :: amold
  class(psb_c_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  ! Local Variables
  integer(psb_ipk_)  :: ictxt, me,np
  integer(psb_ipk_)  :: err,i,k, err_act, iszv, newsz, nplevs, mxplevs
  real(psb_spk_)     :: mnaggratio
  integer(psb_ipk_)  :: coarse_solve_id
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'amg_c_smoothers_bld'
  info = psb_success_
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  !
  if (.not.allocated(prec%precv)) then 
    !! Error: should have called amg_cprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  iszv       = size(prec%precv)
  call psb_bcast(ictxt,iszv)
  if (iszv /= size(prec%precv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if

  if (iszv < 1) then
    ! We should never get here.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif
  
  !
  ! Issue a warning for inconsistent changes to COARSE_SOLVE
  ! but only if it really is a multilevel
  !
  if ((me == psb_root_).and.(iszv>1)) then 
    coarse_solve_id = prec%precv(iszv)%parms%coarse_solve
    select case (coarse_solve_id)
    case(amg_umf_,amg_slu_)
      if (prec%precv(iszv)%sm%sv%get_id() /= coarse_solve_id) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id)
        if (prec%precv(iszv)%parms%coarse_mat == amg_repl_mat_) then
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else if (prec%precv(iszv)%parms%coarse_mat == amg_distr_mat_) then
          write(psb_err_unit,*) ' but I am building BJAC with ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        end if
        write(psb_err_unit,*) 'This may happen if: '
        write(psb_err_unit,*) '  1. coarse_subsolve has been reset, or '
        write(psb_err_unit,*) '  2. the solver ', amg_fact_names(coarse_solve_id),&
             & ' was not configured at MLD2P4 build time, or'
        write(psb_err_unit,*) '  3. an unsupported solver setup was specified.'
      end if
      if (prec%precv(iszv)%parms%coarse_mat /= amg_repl_mat_) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse matrix was requested as replicated', &
             & ' but it has been changed to distributed.'
      end if
            
    case(psb_ilu_n_, psb_ilu_t_,psb_milu_n_)
      if (prec%precv(iszv)%sm%sv%get_id() /= psb_ilu_n_) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id)
        if (prec%precv(iszv)%parms%coarse_mat == amg_repl_mat_) then
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else if (prec%precv(iszv)%parms%coarse_mat == amg_distr_mat_) then
          write(psb_err_unit,*) ' but I am building BJAC with ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        end if
        write(psb_err_unit,*) &
             &'This may happen if coarse_subsolve has been reset'
      end if
      if (prec%precv(iszv)%parms%coarse_mat /= amg_repl_mat_) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id),&
             & ' but the coarse matrix has been changed to distributed'
      end if
      
    case(amg_mumps_)
      if (prec%precv(iszv)%sm%sv%get_id() /= amg_mumps_) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id)
        if (prec%precv(iszv)%parms%coarse_mat == amg_repl_mat_) then
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else if (prec%precv(iszv)%parms%coarse_mat == amg_distr_mat_) then
          write(psb_err_unit,*) ' but I am building BJAC with ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        end if
        write(psb_err_unit,*) 'This may happen if: '
        write(psb_err_unit,*) '  1. coarse_subsolve has been reset, or '
        write(psb_err_unit,*) '  2. the solver ', amg_fact_names(coarse_solve_id),&
             & ' was not configured at MLD2P4 build time, or'
        write(psb_err_unit,*) '  3. an unsupported solver setup was specified.'
      end if
      
    case(amg_sludist_)
      if (prec%precv(iszv)%sm%sv%get_id() /= coarse_solve_id) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id)
        if (prec%precv(iszv)%parms%coarse_mat == amg_repl_mat_) then
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else if (prec%precv(iszv)%parms%coarse_mat == amg_distr_mat_) then
          write(psb_err_unit,*) ' but I am building BJAC with ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        else
          write(psb_err_unit,*) ' but I am building ',&
               & amg_fact_names(prec%precv(iszv)%sm%sv%get_id())
        end if
        write(psb_err_unit,*) 'This may happen if: '
        write(psb_err_unit,*) '  1. coarse_subsolve has been reset, or '
        write(psb_err_unit,*) '  2. the solver ', amg_fact_names(coarse_solve_id), &
             & ' was not configured at MLD2P4 build time, or'
        write(psb_err_unit,*) '  3. an unsupported solver setup was specified.'
      end if
      if (prec%precv(iszv)%parms%coarse_mat /= amg_distr_mat_) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id),&
             & ' but the coarse matrix has been changed to replicated'
      end if
        
    case(amg_bjac_,amg_l1_bjac_,amg_jac_, amg_l1_jac_, amg_gs_, amg_fbgs_, amg_l1_gs_,amg_l1_fbgs_)
      if (prec%precv(iszv)%parms%coarse_mat /= amg_distr_mat_) then
        write(psb_err_unit,*) &
             & 'MLD2P4: Warning: original coarse solver was requested as ',&
             & amg_fact_names(coarse_solve_id),&
             & ' but the coarse matrix has been changed to replicated'
      end if
      
    case default
      ! We should never get here.
      info=psb_err_from_subroutine_
      ch_err='unkn coarse_solve'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
      
    end select
  end if

  ! Sanity check: need to ensure that the MUMPS local/global NZ
  ! are handled correctly; this is controlled by local vs global solver.
  ! From this point of view, REPL is LOCAL because it owns everyting.
  ! Should really find a better way of handling this. 
  if (prec%precv(iszv)%parms%coarse_mat == amg_repl_mat_) &
       &  call prec%precv(iszv)%sm%sv%set('MUMPS_LOC_GLOB', amg_local_solver_,info)
  !
  ! Now do the real build.
  !

  do i=1, iszv
    !
    ! build the base preconditioner at level i
    !
    call prec%precv(i)%bld(info,amold=amold,vmold=vmold,imold=imold,ilv=i)
    
    if (info /= psb_success_) then 
      write(ch_err,'(a,i7)') 'Error @ level',i
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err=ch_err)
      goto 9999
    endif

  end do

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Exiting with',iszv,' levels'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_c_smoothers_bld

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
subroutine amg_d_base_onelev_build(lv,info,amold,vmold,imold,ilv)
  use psb_base_mod
  use amg_d_onelev_mod, amg_protect_name => amg_d_base_onelev_build
  implicit none
  class(amg_d_onelev_type), target, intent(inout) :: lv
  integer(psb_ipk_), intent(out) :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  integer(psb_ipk_), intent(in), optional :: ilv
  ! Local
  integer(psb_ipk_)  :: err,i,k, err_act
  integer(psb_ipk_)  :: ictxt, me, np
  integer(psb_ipk_)  :: debug_level, debug_unit
  character(len=20)  :: name, ch_err

  name = 'amg_onelev_build'
  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  if (.not.associated(lv%base_desc)) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='Unassociated base DESC')
    goto 9999
  end if
  info = psb_success_
  ictxt = lv%base_desc%get_ctxt()
  call psb_info(ictxt,me,np)

  if (.not.allocated(lv%sm)) then 
    !! Error: should have called amg_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (.not.allocated(lv%sm%sv)) then 
    !! Error: should have called amg_dprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if
  lv%ac_nz_loc = lv%ac%get_nzeros()
  lv%ac_nz_tot = lv%ac_nz_loc
  select case(lv%parms%coarse_mat)
  case(amg_distr_mat_) 
    call psb_sum(ictxt,lv%ac_nz_tot)
  case(amg_repl_mat_)
    ! Do nothing
  case default
    ! Should never get here
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Wrong lv%parms')
    goto 9999
  end select
  
  
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Calling mlprcbld at level  ',i
  call amg_check_def(lv%parms%sweeps_pre,&
       & 'Jacobi sweeps',izero,is_int_non_negative)
  call amg_check_def(lv%parms%sweeps_post,&
       & 'Jacobi sweeps',izero,is_int_non_negative)

  call lv%sm%build(lv%base_a,lv%base_desc,info)
  if (info == 0) then
    if (allocated(lv%sm2a)) then 
      call lv%sm2a%build(lv%base_a,lv%base_desc,info)
      lv%sm2 => lv%sm2a
    else
      lv%sm2 => lv%sm
    end if
  end if
  if (info /=0 ) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='Smoother bld error')
    goto 9999
  end if

  if (lv%sm%sv%is_global()) then
    if ((lv%parms%sweeps_pre>1).or.(lv%parms%sweeps_post>1)) then
      lv%parms%sweeps_pre  = 1
      lv%parms%sweeps_post = 1
      if (me == 0) then
        write(debug_unit,*) 
        if (present(ilv)) then
          write(debug_unit,*) 'Warning: the solver "',trim(lv%sm%sv%get_fmt()),&
               & '" at level ',ilv
          write(debug_unit,*) '         is configured as a global solver '
        else
          write(debug_unit,*) 'Warning: the solver "',trim(lv%sm%sv%get_fmt()),&
               & '" is configured as a global solver '
        end if
        write(debug_unit,*) '        Pre and post sweeps at this level reset to 1'
      end if
    end if
  end if
  
  if (any((/present(amold),present(vmold),present(imold)/))) &
       & call lv%cnv(info,amold=amold,vmold=vmold,imold=imold)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_d_base_onelev_build

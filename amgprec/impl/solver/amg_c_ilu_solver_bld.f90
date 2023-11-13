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
subroutine amg_c_ilu_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

  use psb_base_mod
  use amg_c_ilu_solver, amg_protect_name => amg_c_ilu_solver_bld

  Implicit None

  ! Arguments
  type(psb_cspmat_type), intent(in), target           :: a
  Type(psb_desc_type), Intent(inout)                  :: desc_a 
  class(amg_c_ilu_solver_type), intent(inout)         :: sv
  integer(psb_ipk_), intent(out)                      :: info
  type(psb_cspmat_type), intent(in), target, optional :: b
  class(psb_c_base_sparse_mat), intent(in), optional  :: amold
  class(psb_c_base_vect_type), intent(in), optional   :: vmold
  class(psb_i_base_vect_type), intent(in), optional   :: imold
  ! Local variables
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota, psb_fctype
!!$    complex(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act, debug_unit, debug_level
  character(len=20)   :: name='c_ilu_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()

  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()
  if (present(b)) then 
    nztota = nztota + b%get_nzeros()
  end if

  call sv%l%csall(n_row,n_row,info,nztota)
  if (info == psb_success_) call sv%u%csall(n_row,n_row,info,nztota)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (allocated(sv%d)) then 
    if (size(sv%d) < n_row) then 
      deallocate(sv%d)
    endif
  endif
  if (.not.allocated(sv%d))  allocate(sv%d(n_row),stat=info)

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  endif

  select case(sv%fact_type)
  case (amg_ilu_n_)
    psb_fctype = amg_ilu_n_
  case (amg_milu_n_)
    psb_fctype = amg_milu_n_
  case (amg_ilu_t_)
    psb_fctype = amg_ilu_t_
  case default
    ! If we end up here, something was wrong up in the call chain. 
    info = psb_err_input_value_invalid_i_
    call psb_errpush(psb_err_input_value_invalid_i_,name,&
         & i_err=(/ithree,sv%fact_type,izero,izero,izero/))
    goto 9999
  end select
  
  select case(sv%fact_type)

  case (amg_ilu_t_)
    !
    ! ILU(k,t)
    !
    select case(sv%fill_in)

    case(:-1) 
      ! Error: fill-in <= -1
      call psb_errpush(psb_err_input_value_invalid_i_,&
           & name,i_err=(/ithree,sv%fill_in,izero,izero,izero/))
      goto 9999

    case(0:)
      ! Fill-in >= 0
      call psb_ilut_fact(sv%fill_in,sv%thresh,&
           & a, sv%l,sv%u,sv%d,info,blck=b)
    end select
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_ilut_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(amg_ilu_n_,amg_milu_n_) 
    !
    ! ILU(k) and MILU(k)
    !
    select case(sv%fill_in)
    case(:-1) 
      ! Error: fill-in <= -1
      call psb_errpush(psb_err_input_value_invalid_i_,&
           & name,i_err=(/ithree,sv%fill_in,izero,izero,izero/))
      goto 9999
    case(0)
      ! Fill-in 0
      ! Separate implementation of ILU(0) for better performance.
      ! There seems to be a problem with the separate implementation of MILU(0),
      ! contained into psb_ilu0_fact. This must be investigated. For the time being,
      ! resort to the implementation of MILU(k) with k=0.
      if (sv%fact_type == amg_ilu_n_) then 
        call psb_ilu0_fact(psb_fctype,a,sv%l,sv%u,&
             & sv%d,info,blck=b)
      else
        call psb_iluk_fact(sv%fill_in,psb_fctype,&
             & a,sv%l,sv%u,sv%d,info,blck=b)
      endif
    case(1:)
      ! Fill-in >= 1
      ! The same routine implements both ILU(k) and MILU(k)
      call psb_iluk_fact(sv%fill_in,psb_fctype,&
           & a,sv%l,sv%u,sv%d,info,blck=b)
    end select
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_iluk_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case default
    ! If we end up here, something was wrong up in the call chain. 
    info = psb_err_input_value_invalid_i_
    call psb_errpush(psb_err_input_value_invalid_i_,name,&
         & i_err=(/ithree,sv%fact_type,izero,izero,izero/))
    goto 9999

  end select

  call sv%l%set_asb()
  call sv%l%trim()
  call sv%u%set_asb()
  call sv%u%trim()
  call sv%dv%bld(sv%d,mold=vmold)

  if (present(amold)) then 
    call sv%l%cscnv(info,mold=amold)
    call sv%u%cscnv(info,mold=amold)
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_c_ilu_solver_bld

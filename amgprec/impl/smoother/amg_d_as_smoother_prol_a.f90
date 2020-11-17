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
subroutine amg_d_as_smoother_prol_a(sm,x,trans,work,info,data)
  use psb_base_mod
  use amg_d_as_smoother, amg_protect_nam => amg_d_as_smoother_prol_a
  implicit none 
  class(amg_d_as_smoother_type), intent(inout) :: sm
  real(psb_dpk_), intent(inout)                 :: x(:)
  character(len=1),intent(in)                    :: trans
  real(psb_dpk_),target, intent(inout)          :: work(:)
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: data
  !Local
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act,isz,int_err(5), data_
  character           :: trans_
  character(len=20)   :: name='d_as_smther_prol_a', ch_err

  call psb_erractionsave(err_act)

  info  = psb_success_
  ctxt = sm%desc_data%get_context()
  call psb_info(ctxt,me,np)

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T')
  case('C')
  case default
    info = psb_err_iarg_invalid_i_
    call psb_errpush(info,name)
    goto 9999
  end select

  if (present(data)) then
    data_ = data
  else
    data_ = psb_comm_ext_
  end if

  if (.not.allocated(sm%sv)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if



  select case(trans_)
  case('N')

    select case (sm%prol) 

    case(psb_none_)
      ! 
      ! Would work anyway, but since it is supposed to do nothing ...
      !        call psb_ovrl(x,sm%desc_data,info,&
      !             & update=sm%prol,work=work)


    case(psb_sum_,psb_avg_) 
      !
      ! Update the overlap of x
      !
      call psb_ovrl(x,sm%desc_data,info,&
           & update=sm%prol,work=work)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_ovrl'
        goto 9999
      end if

    case default
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Invalid amg_sub_prol_')
      goto 9999
    end select

  case('T','C')
    !
    ! With transpose, we have to do it here
    ! 
    if (sm%restr == psb_halo_) then 
      call psb_ovrl(x,sm%desc_data,info,&
           & update=psb_sum_,work=work)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_ovrl'
        goto 9999
      end if
    else if (sm%restr /= psb_none_) then 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Invalid amg_sub_restr_')
      goto 9999
    end if

  case default
    info=psb_err_iarg_invalid_i_
    int_err(1)=6
    ch_err(2:2)=trans
    goto 9999
  end select


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_d_as_smoother_prol_a



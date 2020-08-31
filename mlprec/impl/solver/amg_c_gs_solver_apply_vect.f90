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
subroutine amg_c_gs_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
     & trans,work,wv,info,init,initu)
  
  use psb_base_mod
  use amg_c_gs_solver, amg_protect_name => amg_c_gs_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)               :: desc_data
  class(amg_c_gs_solver_type), intent(inout) :: sv
  type(psb_c_vect_type),intent(inout)         :: x
  type(psb_c_vect_type),intent(inout)         :: y
  complex(psb_spk_),intent(in)                    :: alpha,beta
  character(len=1),intent(in)                   :: trans
  complex(psb_spk_),target, intent(inout)         :: work(:)
  type(psb_c_vect_type),intent(inout)         :: wv(:)
  integer(psb_ipk_), intent(out)                :: info
  character, intent(in), optional                :: init
  type(psb_c_vect_type),intent(inout), optional   :: initu

  integer(psb_ipk_)   :: n_row,n_col, itx, itxst
  complex(psb_spk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  complex(psb_spk_), allocatable :: temp(:)
  integer(psb_ipk_)   :: ictxt,np,me,i, err_act
  character          :: trans_, init_
  character(len=20)  :: name='c_gs_solver_apply'

  call psb_erractionsave(err_act)
  ictxt = desc_data%get_ctxt()
  call psb_info(ictxt,me,np)
  info = psb_success_

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
!!$  case('T')
!!$  case('C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select
  
  if (present(init)) then
    init_ = psb_toupper(init)
  else
    init_='Z'
  end if


  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()


  if (x%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,&
         & i_err=(/itwo,n_row,izero,izero,izero/))
    goto 9999
  end if
  if (y%get_nrows() < n_row) then 
    info = 36
    call psb_errpush(info,name,& 
         & i_err=(/ithree,n_row,izero,izero,izero/))
    goto 9999
  end if



  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
  endif

  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,&
         & i_err=(/5*n_col,izero,izero,izero,izero/),&
         & a_err='complex(psb_spk_)')
    goto 9999      
  end if

  if (size(wv) < 2) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='invalid wv size')
    goto 9999
  end if

  associate(tw => wv(1), xit => wv(2))
    itxst = 1 
    select case (init_)
    case('Z') 
      call psb_geaxpby(cone,x,czero,tw,desc_data,info)
      call psb_spsm(cone,sv%l,tw,czero,xit,desc_data,info)
      itxst = 2
    case('Y')
      call psb_geaxpby(cone,y,czero,xit,desc_data,info)
    case('U')
      if (.not.present(initu)) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='missing initu to smoother_apply')
        goto 9999
      end if
      call psb_geaxpby(cone,initu,czero,xit,desc_data,info)
    case default
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong  init to smoother_apply')
      goto 9999
    end select

    select case(trans_)
    case('N')
      if (sv%eps <=szero) then
        !
        ! Fixed number of iterations
        !
        !
        do itx=itxst,sv%sweeps
          call psb_geaxpby(cone,x,czero,tw,desc_data,info)
          ! Update with U. The off-diagonal block is taken care
          ! from the Jacobi smoother, hence this is purely local. 
          call psb_spmm(-cone,sv%u,xit,cone,tw,desc_data,info,doswap=.false.)
          call psb_spsm(cone,sv%l,tw,czero,xit,desc_data,info)
        end do

        call psb_geaxpby(alpha,xit,beta,y,desc_data,info)

      else
        !
        ! Iterations to convergence, not implemented right now. 
        !
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='EPS>0 not implemented in GS subsolve')
        goto 9999

      end if

    case default
      info = psb_err_internal_error_
      call psb_errpush(info,name,& 
           & a_err='Invalid TRANS in GS subsolve')
      goto 9999
    end select


    if (info /= psb_success_) then

      call psb_errpush(psb_err_internal_error_,name,& 
           & a_err='Error in subsolve')
      goto 9999
    endif
  end associate
  
  if (n_col <= size(work)) then 
    if ((4*n_col+n_col) <= size(work)) then 
    else
      deallocate(aux)
    endif
  else
    deallocate(ww,aux)
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_c_gs_solver_apply_vect

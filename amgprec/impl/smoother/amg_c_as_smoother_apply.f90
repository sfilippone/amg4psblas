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
subroutine amg_c_as_smoother_apply(alpha,sm,x,beta,y,desc_data,trans,&
     & sweeps,work,info,init,initu)
  use psb_base_mod
  use amg_c_as_smoother, amg_protect_nam => amg_c_as_smoother_apply
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(amg_c_as_smoother_type), intent(inout) :: sm
  complex(psb_spk_),intent(inout)         :: x(:)
  complex(psb_spk_),intent(inout)         :: y(:)
  complex(psb_spk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)           :: trans
  integer(psb_ipk_), intent(in)         :: sweeps
  complex(psb_spk_),target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(out)        :: info
  character, intent(in), optional       :: init
  complex(psb_spk_),intent(inout), optional :: initu(:)

  integer(psb_ipk_)  :: n_row,n_col, nrow_d, i
  complex(psb_spk_), pointer :: aux(:)
  complex(psb_spk_), allocatable  :: tx(:),ty(:), ww(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act,isz,int_err(5)
  character           :: trans_, init_
  character(len=20)   :: name='c_as_smoother_apply', ch_err

  call psb_erractionsave(err_act)

  info  = psb_success_
  ctxt = desc_data%get_context()
  call psb_info(ctxt,me,np)

  if (present(init)) then
    init_ = psb_toupper(init)
  else
    init_='Z'
  end if

  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N')
  case('T')
  case('C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select

  if (.not.allocated(sm%sv)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if


  n_row  = sm%desc_data%get_local_rows()
  n_col  = sm%desc_data%get_local_cols()
  nrow_d = desc_data%get_local_rows()
  isz    = max(n_row,N_COL)

  if ((4*isz) <= size(work)) then 
    aux => work(1:)
  else
    allocate(aux(4*isz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name,&
           & i_err=(/4*isz,izero,izero,izero,izero/),&
           & a_err='complex(psb_spk_)')
      goto 9999      
    end if
  endif

  if ((.not.sm%sv%is_iterative()).and.(sweeps == 1).and.(sm%novr==0)) then 
    !
    ! Shortcut: in this case there is nothing else to be done. 
    !
    call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif

  else if (sweeps >= 0) then 
    !
    !
    ! Apply multiple sweeps of an AS solver
    ! to compute an approximate solution of a linear system.
    !
    !
    call psb_geasb(tx,sm%desc_data,info)
    call psb_geasb(ty,sm%desc_data,info)
    call psb_geasb(ww,sm%desc_data,info)

    !
    !  Unroll  the first iteration and fold it inside SELECT CASE
    !  this will save one SPMM when INIT=Z, and will be
    !  significant when sweeps=1 (a common case)
    !
    call psb_geaxpby(cone,x,czero,tx,desc_data,info)    
    if (info == 0) call sm%apply_restr(tx,trans_,aux,info)
    if (info == 0) call psb_geaxpby(cone,tx,czero,ww,sm%desc_data,info)

    select case (init_)
    case('Z')
      call sm%sv%apply(cone,ww,czero,ty,sm%desc_data,trans_,aux,info,init='Z') 

    case('Y')
      call psb_geaxpby(cone,y,czero,ty,desc_data,info)
      if (info == 0) call sm%apply_restr(ty,trans_,aux,info)
      if (info == 0) call psb_spmm(-cone,sm%nd,ty,cone,ww,sm%desc_data,info,&
           & work=aux,trans=trans_)
      call sm%sv%apply(cone,ww,czero,ty,desc_data,trans_,aux,info,init='Y')             

    case('U')
      if (.not.present(initu)) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='missing initu to smoother_apply')
        goto 9999
      end if
      call psb_geaxpby(cone,initu,czero,ty,desc_data,info)
      if (info == 0) call sm%apply_restr(ty,trans_,aux,info)
      if (info == 0) call psb_spmm(-cone,sm%nd,ty,cone,ww,sm%desc_data,info,&
           & work=aux,trans=trans_)
      call sm%sv%apply(cone,ww,czero,ty,desc_data,trans_,aux,info,init='Y')             

    case default
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong  init to smoother_apply')
      goto 9999
    end select
    if (info == 0) call sm%apply_prol(ty,trans_,aux,info)
    
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif
    
    do i=1, sweeps-1
      !
      ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
      ! block diagonal part and the remaining part of the local matrix
      ! and Y(j) is the approximate solution at sweep j.
      !
      if (info == 0) call psb_geaxpby(cone,tx,czero,ww,sm%desc_data,info)
      if (info == 0) call psb_spmm(-cone,sm%nd,ty,cone,ww,sm%desc_data,info,&
           & work=aux,trans=trans_)
      
      if (info /= psb_success_) exit

      call sm%sv%apply(cone,ww,czero,ty,sm%desc_data,trans_,aux,info,init='Y') 
      
      if (info /= psb_success_) exit
      if (info == 0) call sm%apply_prol(ty,trans_,aux,info)
      
    end do

    if (info /= psb_success_) then 
      info=psb_err_internal_error_
      call psb_errpush(info,name,&
           & a_err='subsolve with Jacobi sweeps > 1')
      goto 9999      
    end if
    
    !
    ! Compute y = beta*y + alpha*ty (ty == K^(-1)*tx)
    !
    call psb_geaxpby(alpha,ty,beta,y,desc_data,info) 
    
    
  else
    
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,&
         & i_err=(/itwo,sweeps,izero,izero,izero/))
    goto 9999

  endif

  
  if (.not.(4*isz <= size(work))) then 
    deallocate(aux,stat=info)
  endif
  if (info ==0) deallocate(ww,tx,ty,stat=info)
  if (info /= 0) then
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_c_as_smoother_apply

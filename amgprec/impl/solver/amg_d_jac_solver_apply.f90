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
subroutine amg_d_jac_solver_apply(alpha,sv,x,beta,y,desc_data,trans,&
     & work,info,init,initu)

  use psb_base_mod
  use amg_d_diag_solver
  use psb_base_krylov_conv_mod, only : log_conv
  use amg_d_jac_solver, amg_protect_name => amg_d_jac_solver_apply
  implicit none
  type(psb_desc_type), intent(in)                 :: desc_data
  class(amg_d_jac_solver_type), intent(inout) :: sv
  real(psb_dpk_),intent(inout)           :: x(:)
  real(psb_dpk_),intent(inout)           :: y(:)
  real(psb_dpk_),intent(in)                      :: alpha,beta
  character(len=1),intent(in)                     :: trans
  real(psb_dpk_),target, intent(inout)           :: work(:)
  integer(psb_ipk_), intent(out)                  :: info
  character, intent(in), optional                :: init
  real(psb_dpk_),intent(inout), optional :: initu(:)
  !
  integer(psb_ipk_)    :: n_row,n_col, sweeps
  real(psb_dpk_), pointer :: aux(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act
  character           :: trans_, init_
  real(psb_dpk_)      :: res, resdenum
  character(len=20)   :: name='d_jac_solver_apply_v'

  call psb_erractionsave(err_act)

  info = psb_success_
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
  case('T','C')
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select


  
  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()
  sweeps = sv%sweeps
  if (4*n_col <= size(work)) then
    aux => work(:)
  else
    allocate(aux(4*n_col),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_request_
      call psb_errpush(info,name,&
           & i_err=(/4*n_col,izero,izero,izero,izero/),&
           & a_err='real(psb_dpk_)')
      goto 9999
    end if
  endif

  if (sweeps >= 0) then
      !
      ! This means we are dealing with a pure Jacobi smoother/solver.
      !
      associate(tx => aux(1:n_col), ty => aux(n_col+1:2*n_col))
        select case (init_)
        case('Z')

          call inner_mlt(n_row,done,sv%dv%v%v,x,dzero,ty,trans=trans_)

        case('Y')
          call psb_geaxpby(done,x,dzero,tx,desc_data,info)
          call psb_geaxpby(done,y,dzero,ty,desc_data,info)
          call psb_spmm(-done,sv%a,ty,done,tx,desc_data,info,&
               & work=aux,trans=trans_, doswap=.false.)
          call inner_mlt(n_row,done,sv%dv%v%v,tx,dzero,ty,trans=trans_)

        case('U')
          if (.not.present(initu)) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='missing initu to smoother_apply')
            goto 9999
          end if
          call psb_geaxpby(done,x,dzero,tx,desc_data,info)
          call psb_geaxpby(done,initu,dzero,ty,desc_data,info)
          call psb_spmm(-done,sv%a,ty,done,tx,desc_data,info,&
               & work=aux,trans=trans_, doswap=.false.)
          call inner_mlt(n_row,done,sv%dv%v%v,tx,dzero,ty,trans=trans_)

        case default
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='wrong  init to smoother_apply')
          goto 9999
        end select

        do i=1, sweeps-1
          !
          ! Compute Y(j+1) =  Y(j)+ D^(-1)*(X-A*Y(j)),
          !  where is the diagonal and  A the matrix.
          !
          call psb_geaxpby(done,x,dzero,tx,desc_data,info)
          call psb_spmm(-done,sv%a,ty,done,tx,desc_data,info,&
               & work=aux,trans=trans_, doswap=.false.)
          if (info /= psb_success_) exit
          call inner_mlt(n_row,done,sv%dv%v%v,tx,done,ty,trans=trans_)
          if (info /= psb_success_) exit
        end do

        if (info == psb_success_) call psb_geaxpby(alpha,ty,beta,y,desc_data,info)

        if (info /= psb_success_) then
          info=psb_err_internal_error_
          call psb_errpush(info,name,&
               & a_err='subsolve with Jacobi sweeps > 1')
          goto 9999
        end if


      end associate


    else
      
      info = psb_err_iarg_neg_
      call psb_errpush(info,name,&
           & i_err=(/itwo,sweeps,izero,izero,izero/))
      goto 9999
      
    end if

  if (.not.(4*n_col <= size(work))) then
    deallocate(aux)
  endif


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
contains
  subroutine inner_mlt(n_row,alpha,d,x,beta,y,trans)
    implicit none
    integer(psb_ipk_),intent(in)   :: n_row
    real(psb_dpk_), intent(inout) :: d(:)
    real(psb_dpk_),intent(inout)  :: x(:)
    real(psb_dpk_),intent(inout)  :: y(:)
    real(psb_dpk_),intent(in)     :: alpha,beta
    character(len=1),intent(in)    :: trans

    integer(psb_ipk_) :: i
    
    if (trans_ == 'C') then 
      if (beta == dzero) then 

        if (alpha == dzero) then 
          y(1:n_row) = dzero
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = (d(i)) * x(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -(d(i)) * x(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * (d(i)) * x(i)
          end do
        end if

      else if (beta == done) then 

        if (alpha == dzero) then 
          !y(1:n_row) = dzero
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = (d(i)) * x(i) + y(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -(d(i)) * x(i)  + y(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * (d(i)) * x(i) + y(i)
          end do
        end if

      else if (beta == -done) then 

        if (alpha == dzero) then 
          y(1:n_row) = -y(1:n_row)        
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = (d(i)) * x(i) - y(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -(d(i)) * x(i)  - y(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * (d(i)) * x(i) - y(i)
          end do
        end if

      else

        if (alpha == dzero) then 
          y(1:n_row) = beta *y(1:n_row)        
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = (d(i)) * x(i) + beta*y(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -(d(i)) * x(i)  + beta*y(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * (d(i)) * x(i) + beta*y(i)
          end do
        end if

      end if

    else if (trans_ /= 'C') then 

      if (beta == dzero) then 

        if (alpha == dzero) then 
          y(1:n_row) = dzero
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = d(i) * x(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -d(i) * x(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * d(i) * x(i)
          end do
        end if

      else if (beta == done) then 

        if (alpha == dzero) then 
          !y(1:n_row) = dzero
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = d(i) * x(i) + y(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -d(i) * x(i)  + y(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * d(i) * x(i) + y(i)
          end do
        end if

      else if (beta == -done) then 

        if (alpha == dzero) then 
          y(1:n_row) = -y(1:n_row)        
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = d(i) * x(i) - y(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -d(i) * x(i)  - y(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * d(i) * x(i) - y(i)
          end do
        end if

      else

        if (alpha == dzero) then 
          y(1:n_row) = beta *y(1:n_row)        
        else if (alpha == done) then 
          do i=1, n_row
            y(i) = d(i) * x(i) + beta*y(i)
          end do
        else if (alpha == -done) then 
          do i=1, n_row
            y(i) = -d(i) * x(i)  + beta*y(i)
          end do
        else
          do i=1, n_row
            y(i) = alpha * d(i) * x(i) + beta*y(i)
          end do
        end if

      end if

    end if

  end subroutine inner_mlt
end subroutine amg_d_jac_solver_apply

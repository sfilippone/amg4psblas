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
subroutine amg_s_jac_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,&
     & sweeps,work,wv,info,init,initu)

  use psb_base_mod
  use amg_s_diag_solver
  use psb_base_krylov_conv_mod, only : log_conv
  use amg_s_jac_smoother, amg_protect_name => amg_s_jac_smoother_apply_vect
  implicit none
  type(psb_desc_type), intent(in)                 :: desc_data
  class(amg_s_jac_smoother_type), intent(inout) :: sm
  type(psb_s_vect_type),intent(inout)           :: x
  type(psb_s_vect_type),intent(inout)           :: y
  real(psb_spk_),intent(in)                      :: alpha,beta
  character(len=1),intent(in)                     :: trans
  integer(psb_ipk_), intent(in)                   :: sweeps
  real(psb_spk_),target, intent(inout)           :: work(:)
  type(psb_s_vect_type),intent(inout)           :: wv(:)
  integer(psb_ipk_), intent(out)                  :: info
  character, intent(in), optional                :: init
  type(psb_s_vect_type),intent(inout), optional   :: initu
  !
  integer(psb_ipk_)    :: n_row,n_col
  type(psb_s_vect_type)  :: tx, ty, r
  real(psb_spk_), pointer :: aux(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act
  character           :: trans_, init_
  real(psb_dpk_)      :: res, resdenum
  character(len=20)   :: name='s_jac_smoother_apply_v'

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

  if (.not.allocated(sm%sv)) then
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (4*n_col <= size(work)) then
    aux => work(:)
  else
    allocate(aux(4*n_col),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_request_
      call psb_errpush(info,name,&
           & i_err=(/4*n_col,izero,izero,izero,izero/),&
           & a_err='real(psb_spk_)')
      goto 9999
    end if
  endif

  if(sm%checkres) then
    call psb_geall(r,desc_data,info)
    call psb_geasb(r,desc_data,info)
    resdenum = psb_genrm2(x,desc_data,info)
  end if

  if ((.not.sm%sv%is_iterative()).and.((sweeps == 1).or.(sm%nd_nnz_tot==0))) then
    !  if .not.sv%is_iterative, there's no need to pass init
    call sm%sv%apply(alpha,x,beta,y,desc_data,trans_,aux,wv,info)

    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,&
           & name,a_err='Error in sub_aply Jacobi Sweeps = 1')
      goto 9999
    endif

  else if (sweeps >= 0) then
    select type (smsv => sm%sv)
    class is (amg_s_diag_solver_type)
      !
      ! This means we are dealing with a pure Jacobi smoother/solver.
      !
      associate(tx => wv(1), ty => wv(2))
        select case (init_)
        case('Z')

          call sm%sv%apply(sone,x,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Z')

        case('Y')
          call psb_geaxpby(sone,x,szero,tx,desc_data,info)
          call psb_geaxpby(sone,y,szero,ty,desc_data,info)
          call psb_spmm(-sone,sm%pa,ty,sone,tx,desc_data,info,work=aux,trans=trans_)
          call sm%sv%apply(sone,tx,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Y')

        case('U')
          if (.not.present(initu)) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='missing initu to smoother_apply')
            goto 9999
          end if
          call psb_geaxpby(sone,x,szero,tx,desc_data,info)
          call psb_geaxpby(sone,initu,szero,ty,desc_data,info)
          call psb_spmm(-sone,sm%pa,ty,sone,tx,desc_data,info,work=aux,trans=trans_)
          call sm%sv%apply(sone,tx,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Y')

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
          call psb_geaxpby(sone,x,szero,tx,desc_data,info)
          call psb_spmm(-sone,sm%pa,ty,sone,tx,desc_data,info,work=aux,trans=trans_)

          if (info /= psb_success_) exit

          call sm%sv%apply(sone,tx,sone,ty,desc_data,trans_,aux,wv(3:),info,init='Y')

          if (info /= psb_success_) exit

          if ( sm%checkres.and.(mod(i,sm%checkiter) == 0) ) then
            call psb_geaxpby(sone,x,szero,r,r,desc_data,info)
            call psb_spmm(-sone,sm%pa,ty,sone,r,desc_data,info)
            res  = psb_genrm2(r,desc_data,info)
            if( sm%printres ) then
              call log_conv("BJAC",me,i,sm%printiter,res,resdenum,sm%tol)
            end if
            if ( res < sm%tol*resdenum ) then
              if( (sm%printres).and.(mod(sm%printiter,sm%checkiter)/=0) ) &
               & call log_conv("BJAC",me,i,1,res,resdenum,sm%tol)
              exit
            end if
          end if

        end do


        if (info == psb_success_) call psb_geaxpby(alpha,ty,beta,y,desc_data,info)

        if (info /= psb_success_) then
          info=psb_err_internal_error_
          call psb_errpush(info,name,&
               & a_err='subsolve with Jacobi sweeps > 1')
          goto 9999
        end if


      end associate

    class default
      !
      !
      ! Apply multiple sweeps of a block-Jacobi solver
      ! to compute an approximate solution of a linear system.
      !
      !
      if (size(wv) < 2) then
        info = psb_err_internal_error_
        call psb_errpush(info,name,&
             & a_err='invalid wv size in smoother_apply')
        goto 9999
      end if
      associate(tx => wv(1), ty => wv(2))

        !
        !  Unroll  the first iteration and fold it inside SELECT CASE
        !  this will save one AXPBY and one SPMM when INIT=Z, and will be
        !  significant when sweeps=1 (a common case)
        !
        select case (init_)
        case('Z')

          call sm%sv%apply(sone,x,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Z')

        case('Y')
          call psb_geaxpby(sone,x,szero,tx,desc_data,info)
          call psb_geaxpby(sone,y,szero,ty,desc_data,info)
          call psb_spmm(-sone,sm%nd,ty,sone,tx,desc_data,info,work=aux,trans=trans_)
          call sm%sv%apply(sone,tx,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Y')

        case('U')
          if (.not.present(initu)) then
            call psb_errpush(psb_err_internal_error_,name,&
                 & a_err='missing initu to smoother_apply')
            goto 9999
          end if
          call psb_geaxpby(sone,x,szero,tx,desc_data,info)
          call psb_geaxpby(sone,initu,szero,ty,desc_data,info)
          call psb_spmm(-sone,sm%nd,ty,sone,tx,desc_data,info,work=aux,trans=trans_)
          call sm%sv%apply(sone,tx,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Y')

        case default
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='wrong  init to smoother_apply')
          goto 9999
        end select

        do i=1, sweeps-1
          !
          ! Compute Y(j+1) = D^(-1)*(X-ND*Y(j)), where D and ND are the
          ! block diagonal part and the remaining part of the local matrix
          ! and Y(j) is the approximate solution at sweep j.
          !
          call psb_geaxpby(sone,x,szero,tx,desc_data,info)
          call psb_spmm(-sone,sm%nd,ty,sone,tx,desc_data,info,work=aux,trans=trans_)

          if (info /= psb_success_) exit

          call sm%sv%apply(sone,tx,szero,ty,desc_data,trans_,aux,wv(3:),info,init='Y')

          if (info /= psb_success_) exit

          if ( sm%checkres.and.(mod(i,sm%checkiter) == 0) ) then
            call psb_geaxpby(sone,x,szero,r,r,desc_data,info)
            call psb_spmm(-sone,sm%pa,ty,sone,r,desc_data,info)
            res  = psb_genrm2(r,desc_data,info)
            if( sm%printres ) then
              call log_conv("BJAC",me,i,sm%printiter,res,resdenum,sm%tol)
            end if
            if (res < sm%tol*resdenum ) then
              if( (sm%printres).and.( mod(sm%printiter,sm%checkiter) /=0 ) ) &
                & call log_conv("BJAC",me,i,1,res,resdenum,sm%tol)
              exit
            end if
          end if

        end do

        if (info == psb_success_) call psb_geaxpby(alpha,ty,beta,y,desc_data,info)

        if (info /= psb_success_) then
          info=psb_err_internal_error_
          call psb_errpush(info,name,&
               & a_err='subsolve with Jacobi sweeps > 1')
          goto 9999
        end if

      end associate
    end select

  else

    info = psb_err_iarg_neg_
    call psb_errpush(info,name,&
         & i_err=(/itwo,sweeps,izero,izero,izero/))
    goto 9999

  endif

  if (.not.(4*n_col <= size(work))) then
    deallocate(aux)
  endif

  if(sm%checkres) then
    call psb_gefree(r,desc_data,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_s_jac_smoother_apply_vect

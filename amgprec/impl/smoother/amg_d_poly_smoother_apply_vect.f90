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
subroutine amg_d_poly_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,&
     & sweeps,work,wv,info,init,initu)

  use psb_base_mod
  use amg_d_diag_solver
  use psb_base_krylov_conv_mod, only : log_conv
  use amg_d_poly_smoother, amg_protect_name => amg_d_poly_smoother_apply_vect
  implicit none
  type(psb_desc_type), intent(in)                 :: desc_data
  class(amg_d_poly_smoother_type), intent(inout) :: sm
  type(psb_d_vect_type),intent(inout)           :: x
  type(psb_d_vect_type),intent(inout)           :: y
  real(psb_dpk_),intent(in)                      :: alpha,beta
  character(len=1),intent(in)                     :: trans
  integer(psb_ipk_), intent(in)                   :: sweeps
  real(psb_dpk_),target, intent(inout)           :: work(:)
  type(psb_d_vect_type),intent(inout)           :: wv(:)
  integer(psb_ipk_), intent(out)                  :: info
  character, intent(in), optional                :: init
  type(psb_d_vect_type),intent(inout), optional   :: initu
  !
  integer(psb_ipk_)    :: n_row,n_col
  type(psb_d_vect_type)  :: tx, ty, tz, r
  real(psb_dpk_), pointer :: aux(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act
  character           :: trans_, init_
  real(psb_dpk_)      :: res, resdenum
  character(len=20)   :: name='d_poly_smoother_apply_v'

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
           & a_err='real(psb_dpk_)')
      goto 9999
    end if
  endif

  if (size(wv) < 4) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='invalid wv size in smoother_apply')
    goto 9999
  end if
  sm%pdegree = sweeps
  associate(tx => wv(1), ty => wv(2), tz => wv(3), r => wv(4))

    call psb_geaxpby(done,x,dzero,r,desc_data,info)
    call tx%zero()
    call ty%zero()
    call tz%zero()       

    select case(sm%variant)
    case(amg_poly_lottes_)
      block 
        real(psb_dpk_)      :: cz, cr
        !  b == x 
        !  x == tx
        ! 
        do i=1, sweeps
          !   B r_{k-1}
          call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
          cz = (2*i*done-3)/(2*i*done+done)
          cr = (8*i*done-4)/((2*i*done+done)*sm%rho_ba)
          ! z_k =  cz z_{k-1} + cr ty = cz z_{k-1} + cr Br_{k-1}
          call psb_geaxpby(cr,ty,cz,tz,desc_data,info)
          ! r_k =    b-Ax_k  = x -A tx
          call psb_geaxpby(done,tz,done,tx,desc_data,info)
          if (.false.) then
            call psb_geaxpby(done,x,dzero,r,desc_data,info)
            call psb_spmm(-done,sm%pa,tx,done,r,desc_data,info,work=aux,trans=trans_)
          else
            call psb_spmm(-done,sm%pa,tz,done,r,desc_data,info,work=aux,trans=trans_)
          end if
!!$          res  = psb_genrm2(r,desc_data,info)
!!$          write(0,*) 'Polynomial smoother LOTTES ',i,res
          ! x_k =  x_{k-1} + z_k
        end do
      end block

    case(amg_poly_lottes_beta_)

      block 
        real(psb_dpk_)      :: cz, cr
        !  b == x 
        !  x == tx
        !
        if (allocated(sm%poly_beta)) then
          if (size(sm%poly_beta) /= sweeps) deallocate(sm%poly_beta)
        end if
        if (.not.allocated(sm%poly_beta)) then
          call psb_realloc(sweeps,sm%poly_beta,info)
          sm%poly_beta(1:sweeps) = amg_d_poly_beta_mat(1:sweeps,sweeps)
        end if

        do i=1, sweeps
          !   B r_{k-1}
          call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
          cz = (2*i*done-3)/(2*i*done+done)
          cr = (8*i*done-4)/((2*i*done+done)*sm%rho_ba)
          ! z_k =  cz z_{k-1} + cr ty = cz z_{k-1} + cr Br_{k-1}
          call psb_geaxpby(cr,ty,cz,tz,desc_data,info)
          ! r_k =    b-Ax_k  = x -A tx
          call psb_geaxpby(sm%poly_beta(i),tz,done,tx,desc_data,info)
          if (.false.) then
            call psb_geaxpby(done,x,dzero,r,desc_data,info)
            call psb_spmm(-done,sm%pa,tx,done,r,desc_data,info,work=aux,trans=trans_)
          else
            call psb_spmm(-done,sm%pa,tz,done,r,desc_data,info,work=aux,trans=trans_)
          end if
!!$          res  = psb_genrm2(r,desc_data,info)
!!$          write(0,*) 'Polynomial smoother LOTTES_BETA ',i,res
          ! x_k =  x_{k-1} + z_k
        end do
      end block
      
    case(amg_poly_new_)
      block 
        real(psb_dpk_)      :: sigma, theta, delta, rho_old, rho
        !  b == x 
        !  x == tx
        !
        sm%cf_a = amg_d_poly_a_vect(sweeps)

        theta = (done+sm%cf_a)/2
        delta = (done-sm%cf_a)/2
        sigma = theta/delta
        rho_old = done/sigma
        call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
        call psb_geaxpby((done/sm%rho_ba),ty,dzero,r,desc_data,info)
        call psb_geaxpby((done/theta),r,dzero,tz,desc_data,info)
        ! tz == d
        do i=1, sweeps
          ! x_{k+1} = x_k + d_k
          call psb_geaxpby(done,tz,done,tx,desc_data,info)
          ! 
          ! r_{k-1} = r_k - (1/rho(BA)) B A d_k
          call psb_spmm(done,sm%pa,tz,dzero,ty,desc_data,info,work=aux,trans=trans_)
          call sm%sv%apply(-done,ty,done,r,desc_data,trans_,aux,wv(5:),info,init='Z')

          !
          ! d_{k+1} = (rho rho_old) d_k + 2(rho/delta) r_{k+1}
          rho = done/(2*sigma - rho_old)
          call psb_geaxpby((2*rho/delta),r,(rho*rho_old),tz,desc_data,info)
!!$          res  = psb_genrm2(r,desc_data,info)
!!$          write(0,*) 'Polynomial smoother NEW ',i,res
          ! x_k =  x_{k-1} + z_k
        end do
      end block


    case default
      info=psb_err_internal_error_
      call psb_errpush(info,name,&
           & a_err='wrong polynomial variant')
      goto 9999
    end select

    if (info == psb_success_) call psb_geaxpby(alpha,tx,beta,y,desc_data,info)

    if (info /= psb_success_) then
      info=psb_err_internal_error_
      call psb_errpush(info,name,&
           & a_err='polynomial smoother')
      goto 9999
    end if
  end associate

  if (.not.(4*n_col <= size(work))) then
    deallocate(aux)
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_d_poly_smoother_apply_vect

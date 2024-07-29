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
  integer(psb_ipk_), intent(in)                   :: sweeps! this is ignored here, the polynomial degree dictates the value
  real(psb_dpk_),target, intent(inout)           :: work(:)
  type(psb_d_vect_type),intent(inout)           :: wv(:)
  integer(psb_ipk_), intent(out)                  :: info
  character, intent(in), optional                :: init
  type(psb_d_vect_type),intent(inout), optional   :: initu
  ! Timers
  logical, parameter  :: do_timings=.true.
  integer(psb_ipk_), save  :: poly_1=-1, poly_2=-1, poly_3=-1
  integer(psb_ipk_), save  :: poly_mv=-1, poly_sv=-1, poly_vect=-1
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
  
  if ((do_timings).and.(poly_1==-1))       &
       & poly_1 = psb_get_timer_idx("POLY: Chebychev4")
  if ((do_timings).and.(poly_2==-1))       &
       & poly_2 = psb_get_timer_idx("POLY: OptChebychev4")
  if ((do_timings).and.(poly_3==-1))       &
       & poly_3 = psb_get_timer_idx("POLY: OptChebychev1")
  if ((do_timings).and.(poly_mv==-1))       &
       & poly_mv = psb_get_timer_idx("POLY: spMV")
  if ((do_timings).and.(poly_vect==-1))       &
       & poly_vect = psb_get_timer_idx("POLY: Vectors")
  if ((do_timings).and.(poly_sv==-1))       &
       & poly_sv = psb_get_timer_idx("POLY: solver")
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

  associate(tx => wv(1), ty => wv(2), tz => wv(3), r => wv(4))

    call psb_geaxpby(done,x,dzero,r,desc_data,info)
    call tx%zero()
    call ty%zero()
    call tz%zero()

    select case(sm%variant)
    case(amg_cheb_4_)
      if (do_timings) call psb_tic(poly_1)
      block
        real(psb_dpk_)      :: cz, cr
        !  b == x
        !  x == tx
        !
        do i=1, sm%pdegree-1
          !   B r_{k-1}
          if (do_timings) call psb_tic(poly_sv)
          call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z') ! ty = M^{-1} r
          if (do_timings) call psb_toc(poly_sv)
          cz = (2*i*done-3)/(2*i*done+done)
          cr = (8*i*done-4)/((2*i*done+done)*sm%rho_ba)
          if (do_timings) call psb_tic(poly_vect)
          call psb_upd_xyz(cr,cz,done,done,ty,tz,tx,desc_data,info) ! zk = cz * zk-1 + cr * rk-1 
          if (do_timings) call psb_toc(poly_vect)
          if (do_timings) call psb_tic(poly_mv)
          call psb_spmm(-done,sm%pa,tz,done,r,desc_data,info,work=aux,trans=trans_)
          if (do_timings) call psb_toc(poly_mv)
        end do
        if (do_timings) call psb_tic(poly_sv)
        call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z') ! ty = M^{-1} r
        if (do_timings) call psb_toc(poly_sv)
        cz = (2*sm%pdegree*done-3)/(2*sm%pdegree*done+done)
        cr = (8*sm%pdegree*done-4)/((2*sm%pdegree*done+done)*sm%rho_ba)
        if (do_timings) call psb_tic(poly_vect)
        call psb_upd_xyz(cr,cz,done,done,ty,tz,tx,desc_data,info)
        if (do_timings) call psb_toc(poly_vect)
      end block
      if (do_timings) call psb_toc(poly_1)

    case(amg_cheb_4_opt_)
      if (do_timings) call psb_tic(poly_2)
      block
        real(psb_dpk_)      :: cz, cr
        !  b == x
        !  x == tx
        !
        if (allocated(sm%poly_beta)) then
          if (size(sm%poly_beta) /= sm%pdegree) deallocate(sm%poly_beta)
        end if
        if (.not.allocated(sm%poly_beta)) then
          call psb_realloc(sm%pdegree,sm%poly_beta,info)
          sm%poly_beta(1:sm%pdegree) = amg_d_poly_beta_mat(1:sm%pdegree,sm%pdegree)
        end if

        do i=1, sm%pdegree-1
          !   B r_{k-1}
          if (do_timings) call psb_tic(poly_sv)
          call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
          if (do_timings) call psb_toc(poly_sv)
          cz = (2*i*done-3)/(2*i*done+done)
          cr = (8*i*done-4)/((2*i*done+done)*sm%rho_ba)
          if (do_timings) call psb_tic(poly_vect)
          call psb_upd_xyz(cr,cz,sm%poly_beta(i),done,ty,tz,tx,desc_data,info)
          if (do_timings) call psb_toc(poly_vect)
          if (do_timings) call psb_tic(poly_mv)
          call psb_spmm(-done,sm%pa,tz,done,r,desc_data,info,work=aux,trans=trans_)
          if (do_timings) call psb_toc(poly_mv)
        end do
        call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
        cz = (2*sm%pdegree*done-3)/(2*sm%pdegree*done+done)
        cr = (8*sm%pdegree*done-4)/((2*sm%pdegree*done+done)*sm%rho_ba)
        if (do_timings) call psb_tic(poly_vect)
        call psb_upd_xyz(cr,cz,sm%poly_beta(sm%pdegree),done,ty,tz,tx,desc_data,info)
        if (do_timings) call psb_toc(poly_vect)
      end block
      if (do_timings) call psb_toc(poly_2)
    case(amg_cheb_1_opt_)
      if (do_timings) call psb_tic(poly_3)
      block
        real(psb_dpk_)      :: sigma, theta, delta, rho_old, rho
        !  b == x
        !  x == tx
        !

        theta = (done+sm%cf_a)/2
        delta = (done-sm%cf_a)/2
        sigma = theta/delta
        rho_old = done/sigma
        if (do_timings) call psb_tic(poly_sv)
        call sm%sv%apply(done,r,dzero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
        if (do_timings) call psb_toc(poly_sv)
        call psb_geaxpby((done/sm%rho_ba),ty,dzero,r,desc_data,info)
        if (do_timings) call psb_tic(poly_vect)
        call psb_upd_xyz((done/theta),dzero,done,done,r,tz,tx,desc_data,info)
        if (do_timings) call psb_toc(poly_vect)

        ! tz == d
        do i=1, sm%pdegree-1
          !
          !
          ! r_{k-1} = r_k - (1/rho(BA)) B A d_k
          if (do_timings) call psb_tic(poly_mv)
          call psb_spmm(done,sm%pa,tz,dzero,ty,desc_data,info,work=aux,trans=trans_)
          if (do_timings) call psb_toc(poly_mv)
          if (do_timings) call psb_tic(poly_sv)
          call sm%sv%apply(-(done/sm%rho_ba),ty,done,r,desc_data,trans_,aux,wv(5:),info,init='Z')
          if (do_timings) call psb_toc(poly_sv)
          !
          ! d_{k+1} = (rho rho_old) d_k + 2(rho/delta) r_{k+1}
          rho = done/(2*sigma - rho_old)
          if (do_timings) call psb_tic(poly_vect)
          call psb_upd_xyz((2*rho/delta),(rho*rho_old),done,done,r,tz,tx,desc_data,info)
          if (do_timings) call psb_toc(poly_vect)
          rho_old = rho
        end do
      end block
      if (do_timings) call psb_toc(poly_3)
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

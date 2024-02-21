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
subroutine amg_s_poly_smoother_apply_vect(alpha,sm,x,beta,y,desc_data,trans,&
     & sweeps,work,wv,info,init,initu)

  use psb_base_mod
  use amg_s_diag_solver
  use psb_base_krylov_conv_mod, only : log_conv
  use amg_s_poly_smoother, amg_protect_name => amg_s_poly_smoother_apply_vect
  implicit none
  type(psb_desc_type), intent(in)                 :: desc_data
  class(amg_s_poly_smoother_type), intent(inout) :: sm
  type(psb_s_vect_type),intent(inout)           :: x
  type(psb_s_vect_type),intent(inout)           :: y
  real(psb_spk_),intent(in)                      :: alpha,beta
  character(len=1),intent(in)                     :: trans
  integer(psb_ipk_), intent(in)                   :: sweeps! this is ignored here, the polynomial degree dictates the value
  real(psb_spk_),target, intent(inout)           :: work(:)
  type(psb_s_vect_type),intent(inout)           :: wv(:)
  integer(psb_ipk_), intent(out)                  :: info
  character, intent(in), optional                :: init
  type(psb_s_vect_type),intent(inout), optional   :: initu
  !
  integer(psb_ipk_)    :: n_row,n_col
  type(psb_s_vect_type)  :: tx, ty, tz, r
  real(psb_spk_), pointer :: aux(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, i, err_act
  character           :: trans_, init_
  real(psb_spk_)      :: res, resdenum
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
           & a_err='real(psb_spk_)')
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

    call psb_geaxpby(sone,x,szero,r,desc_data,info)
    call tx%zero()
    call ty%zero()
    call tz%zero()

    select case(sm%variant)
    case(amg_poly_lottes_)
      block
        real(psb_spk_)      :: cz, cr
        !  b == x
        !  x == tx
        !
        do i=1, sm%pdegree
          !   B r_{k-1}
          call sm%sv%apply(sone,r,szero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
          cz = (2*i*sone-3)/(2*i*sone+sone)
          cr = (8*i*sone-4)/((2*i*sone+sone)*sm%rho_ba)
          if (.false.) then 
            ! z_k =  cz z_{k-1} + cr ty = cz z_{k-1} + cr Br_{k-1}
            call psb_geaxpby(cr,ty,cz,tz,desc_data,info)
            ! r_k =    b-Ax_k  = x -A tx
            call psb_geaxpby(sone,tz,sone,tx,desc_data,info)
          else
            call psb_abgdxyz(cr,cz,sone,sone,ty,tz,tx,desc_data,info)
          end if
          if (.false.) then
            call psb_geaxpby(sone,x,szero,r,desc_data,info)
            call psb_spmm(-sone,sm%pa,tx,sone,r,desc_data,info,work=aux,trans=trans_)
          else
            call psb_spmm(-sone,sm%pa,tz,sone,r,desc_data,info,work=aux,trans=trans_)
          end if
!!$        res  = psb_genrm2(r,desc_data,info)
!!$        write(0,*) 'Polynomial smoother  LOTTES',i,res
          ! x_k =  x_{k-1} + z_k
        end do
      end block

    case(amg_poly_lottes_beta_)

      block
        real(psb_spk_)      :: cz, cr
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

        do i=1, sm%pdegree
          !   B r_{k-1}
          call sm%sv%apply(sone,r,szero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
          cz = (2*i*sone-3)/(2*i*sone+sone)
          cr = (8*i*sone-4)/((2*i*sone+sone)*sm%rho_ba)
          if (.false.) then 
            ! z_k =  cz z_{k-1} + cr ty = cz z_{k-1} + cr Br_{k-1}
            call psb_geaxpby(cr,ty,cz,tz,desc_data,info)
            ! r_k =    b-Ax_k  = x -A tx
            call psb_geaxpby(sm%poly_beta(i),tz,sone,tx,desc_data,info)
          else
            call psb_abgdxyz(cr,cz,sm%poly_beta(i),sone,ty,tz,tx,desc_data,info)
          end if
          if (.false.) then
            call psb_geaxpby(sone,x,szero,r,desc_data,info)
            call psb_spmm(-sone,sm%pa,tx,sone,r,desc_data,info,work=aux,trans=trans_)
          else
            call psb_spmm(-sone,sm%pa,tz,sone,r,desc_data,info,work=aux,trans=trans_)
          end if
!!$        res  = psb_genrm2(r,desc_data,info)
!!$        write(0,*) 'Polynomial smoother LOTTES_BETA',i,res
          ! x_k =  x_{k-1} + z_k
        end do
      end block

    case(amg_poly_new_)
      block
        real(psb_spk_)      :: sigma, theta, delta, rho_old, rho
        !  b == x
        !  x == tx
        !
        sm%cf_a = amg_d_poly_a_vect(sm%pdegree)

        theta = (sone+sm%cf_a)/2
        delta = (sone-sm%cf_a)/2
        sigma = theta/delta
        rho_old = sone/sigma
        call sm%sv%apply(sone,r,szero,ty,desc_data,trans_,aux,wv(5:),info,init='Z')
        call psb_geaxpby((sone/sm%rho_ba),ty,szero,r,desc_data,info)
        if (.false.) then 
          call psb_geaxpby((sone/theta),r,szero,tz,desc_data,info)
          call psb_geaxpby(sone,tz,sone,tx,desc_data,info)
        else
          call psb_abgdxyz((sone/theta),szero,sone,sone,r,tz,tx,desc_data,info)
        end if
          
        ! tz == d
        do i=1, sm%pdegree
          !
           !
          ! r_{k-1} = r_k - (1/rho(BA)) B A d_k
          call psb_spmm(sone,sm%pa,tz,szero,ty,desc_data,info,work=aux,trans=trans_)
          call sm%sv%apply(-sone,ty,sone,r,desc_data,trans_,aux,wv(5:),info,init='Z')

          !
          ! d_{k+1} = (rho rho_old) d_k + 2(rho/delta) r_{k+1}
          rho = sone/(2*sigma - rho_old)
          if (.false.) then 
            call psb_geaxpby((2*rho/delta),r,(rho*rho_old),tz,desc_data,info)
            call psb_geaxpby(sone,tz,sone,tx,desc_data,info)
          else
            call psb_abgdxyz((2*rho/delta),(rho*rho_old),sone,sone,r,tz,tx,desc_data,info)
          end if
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

end subroutine amg_s_poly_smoother_apply_vect

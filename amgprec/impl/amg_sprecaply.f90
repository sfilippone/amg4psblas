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
! File: amg_sprecaply.f90
!
! Subroutine: amg_sprecaply
! Version:    real
!
!  This routine applies the preconditioner built by amg_sprecbld, i.e. it computes
!
!                             Y = op(M^(-1)) * X,
!  where
!  - M is the preconditioner,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors.
!  This operation is performed at each iteration of a preconditioned Krylov solver.
!
!
! Arguments:
!    prec       -  type(amg_sprec_type), input.
!                  The preconditioner data structure containing the local part
!                  of the preconditioner to be applied.
!    x          -  real(psb_spk_), dimension(:), input.
!                  The local part of the vector X in Y=op(M^(-1))*X.
!    y          -  real(psb_spk_), dimension(:), output.
!                  The local part of the vector Y in Y=op(M^(-1))*X.
!    desc_data  -  type(psb_desc_type), input.
!                  The communication descriptor associated to the matrix to be
!                  preconditioned.
!    info       -  integer, output.
!                  Error code.
!    trans      -  character(len=1), optional.
!                  If trans='N','n' then op(M^(-1)) = M^(-1);
!                  if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!    work       -  real(psb_spk_), dimension (:), optional, target.
!                  Workspace. Its size must be at
!                  least 4*desc_data%get_local_cols().
!    
subroutine amg_sprecaply(prec,x,y,desc_data,info,trans,work)

  use psb_base_mod
  use amg_s_inner_mod!, amg_protect_name => amg_sprecaply
  
  implicit none
  
  ! Arguments
  type(psb_desc_type),intent(in)    :: desc_data
  type(amg_sprec_type), intent(inout)  :: prec
  real(psb_spk_),intent(inout)   :: x(:)
  real(psb_spk_),intent(inout)   :: y(:)
  integer(psb_ipk_), intent(out)    :: info
  character(len=1), optional        :: trans
  real(psb_spk_),intent(inout), optional, target  :: work(:)

  ! Local variables
  character     :: trans_ 
  real(psb_spk_), pointer :: work_(:)
  real(psb_spk_), allocatable :: w1(:), w2(:) 
  
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_ipk_)   :: err_act,iwsz, k, nswps
  character(len=20)   :: name

  name='amg_sprecaply'
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

  if (present(trans)) then 
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    iwsz = max(1,4*desc_data%get_local_cols())
    allocate(work_(iwsz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name, &
           & i_err=(/iwsz,izero,izero,izero,izero/),&
           & a_err='real(psb_spk_)')
      goto 9999      
    end if

  end if

  if (.not.(allocated(prec%precv))) then 
    !! Error 1: should call amg_sprecbld
    info=3112
    call psb_errpush(info,name)
    goto 9999
  end if

  if (size(prec%precv) >1) then
    !
    ! Number of levels > 1: apply the multilevel preconditioner
    ! 
    call amg_mlprec_aply(sone,prec,x,szero,y,desc_data,trans_,work_,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smlprec_aply')
      goto 9999
    end if

  else  if (size(prec%precv) == 1) then
    !
    ! Number of levels = 1: apply the base preconditioner
    !
    if (allocated(prec%precv(1)%sm2a)) then
      nswps = max(prec%precv(1)%parms%sweeps_pre,prec%precv(1)%parms%sweeps_post)
      !
      ! This is a kludge for handling the symmetrized GS case.
      ! Will need some rethinking. 
      ! 
      call psb_geasb(w1,desc_data,info,scratch=.true.)
      call psb_geasb(w2,desc_data,info,scratch=.true.)

      call psb_geaxpby(sone,x,szero,w1,desc_data,info)
      select case(trans_)
      case ('N')
        do k=1, nswps
          call prec%precv(1)%sm%apply(sone,w1,szero,w2,desc_data,trans_,&
               & ione, work_,info)
          call prec%precv(1)%sm2a%apply(sone,w2,szero,w1,desc_data,trans_,&
               & ione, work_,info)
        end do

      case('T','C')
        do k=1, nswps
          call prec%precv(1)%sm2a%apply(sone,w1,szero,w2,desc_data,trans_,&
               & ione, work_,info)
          call prec%precv(1)%sm%apply(sone,w2,szero,w1,desc_data,trans_,&
               & ione, work_,info)
        end do
      case default
        info = psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='Invalid trans')
        goto 9999         
      end select
      call psb_geaxpby(sone,w1,szero,y,desc_data,info)
      call psb_gefree(w1,desc_data,info)        
      call psb_gefree(w2,desc_data,info)

    else
      nswps = prec%precv(1)%parms%sweeps_pre
      call prec%precv(1)%sm%apply(sone,x,szero,y,desc_data,trans_,&
           & nswps, work_,info)
    end if
  else 
    info = psb_err_from_subroutine_ai_
    call psb_errpush(info,name,a_err='Invalid size of precv',&
         & i_Err=(/ione*size(prec%precv),izero,izero,izero,izero/))
    goto 9999
  endif

  ! If the original distribution has an overlap we should fix that. 
  call psb_halo(y,desc_data,info,data=psb_comm_mov_)


  if (present(work)) then 
  else
    deallocate(work_)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_sprecaply


!
! Subroutine: amg_sprecaply1
! Version:    real
!
!  Applies the preconditioner built by amg_sprecbld, i.e. computes
!
!                             X = op(M^(-1)) * X,
!  where
!  - M is the preconditioner,
!  - op(M^(-1)) is M^(-1) or its transpose, according to the value of trans,
!  - X is a vectors.
!  This operation is performed at each iteration of a preconditioned Krylov solver.
!
!  This routine differs from amg_sprecaply because the preconditioned vector X
!  overwrites the original one.
!
!
! Arguments:
!    prec       -  type(amg_sprec_type), input.
!                  The preconditioner data structure containing the local part
!                  of the preconditioner to be applied.
!    x          -  real(psb_spk_), dimension(:), input/output.
!                  The local part of vector X in X := op(M^(-1)) * X.
!    desc_data  -  type(psb_desc_type), input.
!                  The communication descriptor associated to the matrix to be
!                  preconditioned.
!    info       -  integer, output.
!                  Error code.
!    trans      -  character(len=1), optional.
!                  If trans='N','n' then op(M^(-1)) = M^(-1);
!                  if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!  
subroutine amg_sprecaply1(prec,x,desc_data,info,trans)

  use psb_base_mod
  use amg_s_inner_mod!, amg_protect_name => amg_sprecaply1

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)    :: desc_data
  type(amg_sprec_type), intent(inout)  :: prec
  real(psb_spk_),intent(inout)   :: x(:)
  integer(psb_ipk_), intent(out)    :: info
  character(len=1), optional        :: trans

  ! Local variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_ipk_)   :: err_act
  real(psb_spk_), pointer :: ww(:), w1(:)
  character(len=20)   :: name

  name='amg_sprecaply1'
  info = psb_success_
  call psb_erractionsave(err_act)
  

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

  allocate(ww(size(x)),w1(size(x)),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name, &
         & i_err=(/itwo*size(x),izero,izero,izero,izero/),&
         & a_err='real(psb_spk_)')
    goto 9999      
  end if

  call prec%apply(x,ww,desc_data,info,trans=trans,work=w1)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_precaply')
    goto 9999
  end if

  x(:) = ww(:)
  deallocate(ww,w1,stat=info)
  if (info /= psb_success_) then
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine amg_sprecaply1



subroutine amg_sprecaply2_vect(prec,x,y,desc_data,info,trans,work)

  use psb_base_mod
  use amg_s_inner_mod!, amg_protect_name => amg_sprecaply2_vect
  
  implicit none
  
  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(amg_sprec_type), intent(inout) :: prec
  type(psb_s_vect_type),intent(inout) :: x
  type(psb_s_vect_type),intent(inout) :: y
  integer(psb_ipk_), intent(out)      :: info
  character(len=1), optional          :: trans
  real(psb_spk_),intent(inout), optional, target  :: work(:)

  ! Local variables
  character     :: trans_ 
  real(psb_spk_), pointer :: work_(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_ipk_)   :: err_act,iwsz, k, nswps
  logical             :: do_alloc_wrk
  character(len=20)   :: name

  name='amg_sprecaply'
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

  if (present(trans)) then 
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    iwsz = max(1,4*desc_data%get_local_cols())
    allocate(work_(iwsz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name, &
           & i_err=(/iwsz,izero,izero,izero,izero/),&
           & a_err='real(psb_spk_)')
      goto 9999      
    end if

  end if

  if (.not.(allocated(prec%precv))) then 
    !! Error 1: should call amg_sprecbld
    info=3112
    call psb_errpush(info,name)
    goto 9999
  end if
  
  do_alloc_wrk = .not.allocated(prec%precv(1)%wrk)
  if (do_alloc_wrk) call prec%allocate_wrk(info,vmold=x%v)

  if (size(prec%precv) >1) then
    !
    ! Number of levels > 1: apply the multilevel preconditioner
    !
    ! FIXME: generic name causes an ICE with Intel
    call amg_smlprec_aply_vect(sone,prec,x,szero,y,desc_data,trans_,work_,info)

    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smlprec_aply')
      goto 9999
    end if

  else  if (size(prec%precv) == 1) then
    !
    ! Number of levels = 1: apply the base preconditioner
    !
    nswps = max(prec%precv(1)%parms%sweeps_pre,prec%precv(1)%parms%sweeps_post)

    associate(w1 => prec%precv(1)%wrk%vx2l, w2 => prec%precv(1)%wrk%vy2l,&
         & wv =>  prec%precv(1)%wrk%wv)
      if (allocated(prec%precv(1)%sm2a)) then
        !
        ! This is a kludge for handling the symmetrized GS case.
        ! Will need some rethinking. 
        !
        call psb_geaxpby(sone,x,szero,w1,desc_data,info)
        select case(trans_)
        case ('N')
          do k=1, nswps
            if (info == 0) call prec%precv(1)%sm%apply(sone,w1,szero,w2,desc_data,trans_,&
                 & ione, work_,wv,info)
            if (info == 0) call prec%precv(1)%sm2a%apply(sone,w2,szero,w1,desc_data,trans_,&
                 & ione, work_,wv,info)
          end do
          
        case('T','C')
          do k=1, nswps
            if (info == 0) call prec%precv(1)%sm2a%apply(sone,w1,szero,w2,desc_data,trans_,&
                 & ione, work_,wv,info)
            if (info == 0) call prec%precv(1)%sm%apply(sone,w2,szero,w1,desc_data,trans_,&
                 & ione, work_,wv,info)
          end do
        case default
          info = psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='Invalid trans')
          goto 9999         
        end select
        if (info == 0) call psb_geaxpby(sone,w1,szero,y,desc_data,info)
      else
        if (info == 0) call prec%precv(1)%sm%apply(sone,x,szero,y,desc_data,trans_,&
             & nswps,work_,wv,info)
      end if
    end associate
    if (psb_errstatus_fatal())   info = psb_err_internal_error_
    if (info /= 0) then
      info = psb_err_from_subroutine_ai_
      call psb_errpush(info,name,a_err='Smoother application',&
           & i_Err=(/ione*size(prec%precv),izero,izero,izero,izero/))
      goto 9999
    end if
    
  else 
    
    info = psb_err_from_subroutine_ai_
    call psb_errpush(info,name,a_err='Invalid size of precv',&
         & i_Err=(/ione*size(prec%precv),izero,izero,izero,izero/))
    goto 9999
  endif

  ! If the original distribution has an overlap we should fix that. 
  call psb_halo(y,desc_data,info,data=psb_comm_mov_)

  if (do_alloc_wrk) call prec%free_wrk(info)

  if (present(work)) then 
  else
    deallocate(work_)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_sprecaply2_vect


subroutine amg_sprecaply1_vect(prec,x,desc_data,info,trans,work)

  use psb_base_mod
  use amg_s_inner_mod!, amg_protect_name => amg_sprecaply1_vect

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)      :: desc_data
  type(amg_sprec_type), intent(inout) :: prec
  type(psb_s_vect_type),intent(inout) :: x
  integer(psb_ipk_), intent(out)      :: info
  character(len=1), optional          :: trans
  real(psb_spk_),intent(inout), optional, target  :: work(:)

  ! Local variables
  character     :: trans_ 
  real(psb_spk_), pointer :: work_(:)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_ipk_)   :: err_act,iwsz, k, nswps
  logical             :: do_alloc_wrk
  character(len=20)   :: name

  name='amg_sprecaply'
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

  if (present(trans)) then 
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    iwsz = max(1,4*desc_data%get_local_cols())
    allocate(work_(iwsz),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name, &
           & i_err=(/iwsz,izero,izero,izero,izero/),&
           & a_err='real(psb_spk_)')
      goto 9999      
    end if

  end if

  if (.not.(allocated(prec%precv))) then 
    !! Error 1: should call amg_sprecbld
    info=3112
    call psb_errpush(info,name)
    goto 9999
  end if
  
  do_alloc_wrk = .not.allocated(prec%precv(1)%wrk)
  if (do_alloc_wrk) call prec%allocate_wrk(info,vmold=x%v)

  associate(ww => prec%precv(1)%wrk%vtx, wv =>  prec%precv(1)%wrk%wv)

    if (size(prec%precv) >1) then
      !
      ! Number of levels > 1: apply the multilevel preconditioner
      !
      ! FIXME: generic name causes an ICE with Intel
      call amg_smlprec_aply_vect(sone,prec,x,szero,ww,desc_data,trans_,work_,info)
      if (info == 0) call psb_geaxpby(sone,ww,szero,x,desc_data,info)
      if(info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smlprec_aply')
        goto 9999
      end if

    else  if (size(prec%precv) == 1) then
      !
      ! Number of levels = 1: apply the base preconditioner
      !
      nswps = max(prec%precv(1)%parms%sweeps_pre,prec%precv(1)%parms%sweeps_post)
      if (allocated(prec%precv(1)%sm2a)) then
        !
        ! This is a kludge for handling the symmetrized GS case.
        ! Will need some rethinking. 
        ! 
        select case(trans_)
        case ('N')
          do k=1, nswps
            if (info == 0) call prec%precv(1)%sm%apply(sone,x,szero,ww,desc_data,trans_,&
                 & ione, work_,wv,info)
            if (info == 0) call prec%precv(1)%sm2a%apply(sone,ww,szero,x,desc_data,trans_,&
                 & ione, work_,wv,info)
          end do
        case('T','C')
          do k=1, nswps
            if (info == 0) call prec%precv(1)%sm2a%apply(sone,x,szero,ww,desc_data,trans_,&
                 & ione, work_,wv,info)
            if (info == 0) call prec%precv(1)%sm%apply(sone,ww,szero,x,desc_data,trans_,&
                 & ione, work_,wv,info)
          end do
        case default
          info = psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='Invalid trans')
          goto 9999         
        end select

      else
        if (info == 0) call prec%precv(1)%sm%apply(sone,x,szero,ww,desc_data,trans_,&
             & nswps, work_,wv,info)
        if (info == 0) call psb_geaxpby(sone,ww,szero,x,desc_data,info)
      end if

      if (psb_errstatus_fatal())   info = psb_err_internal_error_
      if (info /=0) then
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='Smoother application',&
             & i_Err=(/ione*size(prec%precv),izero,izero,izero,izero/))
        goto 9999
      end if

    else 

      info = psb_err_from_subroutine_ai_
      call psb_errpush(info,name,a_err='Invalid size of precv',&
           & i_Err=(/ione*size(prec%precv),izero,izero,izero,izero/))
      goto 9999
    endif
  end associate

  ! If the original distribution has an overlap we should fix that. 
  call psb_halo(x,desc_data,info,data=psb_comm_mov_)

  if (do_alloc_wrk) call prec%free_wrk(info)

  if (present(work)) then 
  else
    deallocate(work_)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_sprecaply1_vect

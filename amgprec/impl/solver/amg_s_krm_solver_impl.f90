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
! moved here from amg4psblas-extension
!
!
!                             AMG4PSBLAS  Extensions
!
!    (C) Copyright 2019
!
!                        Salvatore Filippone  Cranfield University
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
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
! File: amg_s_krm_solver_impl.f90
!
!  This is the implementation file corresponding to amg_s_krm_solver_mod.
!
!
subroutine amg_s_krm_solver_bld(a,desc_a,sv,info,b,amold,vmold)

  use psb_base_mod
  use amg_s_krm_solver, amg_protect_name => amg_s_krm_solver_bld

  Implicit None

  ! Arguments
  type(psb_sspmat_type), intent(inout), target        :: a
  Type(psb_desc_type), Intent(inout)                  :: desc_a
  class(amg_s_krm_solver_type), intent(inout)         :: sv
  integer(psb_ipk_), intent(out)                      :: info
  type(psb_sspmat_type), intent(in), target, optional :: b
  class(psb_s_base_sparse_mat), intent(in), optional  :: amold
  class(psb_s_base_vect_type), intent(in), optional   :: vmold
  ! Local variables
  integer(psb_ipk_) :: n_row,n_col, nrow_a, nztota
  integer(psb_lpk_) :: lnr
  integer(psb_ipk_) :: np,me,i, err_act, debug_unit, debug_level
  type(psb_ctxt_type) :: ctxt, l_ctxt
  character(len=20)   :: name='@Z@_krm_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  call sv%prec%free(info)

  if (sv%global) then
    call sv%prec%init(ctxt,sv%kprec,info)
    if (sv%i_sub_solve>0) then
      call sv%prec%set('sub_solve',sv%i_sub_solve,info)
    else
      call sv%prec%set('sub_solve',sv%sub_solve,info)
    end if
    call sv%prec%set('sub_fillin',sv%fillin,info)
    call sv%prec%hierarchy_build(a,desc_a,info)
    call sv%prec%smoothers_build(a,desc_a,info,amold=amold,vmold=vmold)
    sv%a => a
  else
    call psb_init(l_ctxt,np=1_psb_ipk_,basectxt=ctxt,ids=(/me/))
    n_row = desc_a%get_local_rows()
    lnr = n_row
    call psb_cdall(l_ctxt,sv%desc_local,info,mg=lnr,repl=.true.)
    call psb_cdasb(sv%desc_local,info)
    call sv%prec%init(l_ctxt,sv%kprec,info)
    !This is a gigantic kludge, to be revised
    if (sv%i_sub_solve>0) then
      call sv%prec%set('sub_solve',sv%i_sub_solve,info)
    else
      call sv%prec%set('sub_solve',sv%sub_solve,info)
    end if
    call sv%prec%set('sub_fillin',sv%fillin,info)
    call a%csclip(sv%a_local,info,jmax=n_row)
    call sv%prec%hierarchy_build(sv%a_local,sv%desc_local,info)
    call sv%prec%smoothers_build(sv%a_local,sv%desc_local,info,amold=amold,vmold=vmold)
    call psb_geall(sv%x_local,sv%desc_local,info)
    call sv%x_local%zero()
    call psb_geasb(sv%x_local,sv%desc_local,info,mold=vmold)
    call psb_geall(sv%z_local,sv%desc_local,info)
    call sv%z_local%zero()
    call psb_geasb(sv%z_local,sv%desc_local,info,mold=vmold)
  end if

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='inner precbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_s_krm_solver_bld

subroutine amg_s_krm_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
     & trans,work,wv,info,init,initu)

  use psb_base_mod
  use psb_krylov_mod
  use amg_s_krm_solver, amg_protect_name => amg_s_krm_solver_apply_vect

  Implicit None
  type(psb_desc_type), intent(in)             :: desc_data
  class(amg_s_krm_solver_type), intent(inout) :: sv
  type(psb_s_vect_type),intent(inout)         :: x
  type(psb_s_vect_type),intent(inout)         :: y
  real(psb_spk_),intent(in)                    :: alpha,beta
  character(len=1),intent(in)                   :: trans
  real(psb_spk_),target, intent(inout)         :: work(:)
  type(psb_s_vect_type),intent(inout)         :: wv(:)
  integer(psb_ipk_), intent(out)                :: info
  character, intent(in), optional                :: init
  type(psb_s_vect_type),intent(inout), optional   :: initu

  type(psb_s_vect_type)   :: z
  integer(psb_ipk_) :: np,me,i, err_act, debug_unit, debug_level
  type(psb_ctxt_type) :: ctxt
  character(len=20) :: name='@Z@_krm_solver_apply_v', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)

  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_data%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'
  if (sv%global) then
    call psb_geall(z,desc_data,info)
    call z%zero()
    call psb_geasb(z,desc_data,info,mold=x%v)

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' norms ',&
         & psb_genrm2(x,desc_data,info),psb_genrm2(y,desc_data,info),&
         & psb_genrm2(z,desc_data,info)

    call psb_krylov(sv%method,sv%a,sv%prec,x,z,sv%eps,&
         & desc_data,info,itmax=sv%itmax,itrace=sv%itrace,&
         & istop=sv%istopc,irst=sv%irst)
!!$  call sv%prec%apply(x,z,desc_data,info,trans=trans,work=work)
    call psb_geaxpby(alpha,z,beta,y,desc_data,info)
  else
    call psb_geaxpby(sone,x,szero,sv%x_local,sv%desc_local,info)
    call psb_krylov(sv%method,sv%a_local,sv%prec,sv%x_local,sv%z_local,sv%eps,&
         & sv%desc_local,info,itmax=sv%itmax,itrace=sv%itrace,&
         & istop=sv%istopc,irst=sv%irst)
    call psb_geaxpby(alpha,sv%z_local,beta,y,sv%desc_local,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_s_krm_solver_apply_vect


subroutine amg_s_krm_solver_apply(alpha,sv,x,beta,y,desc_data,&
     & trans,work,info,init,initu)
  use psb_base_mod
  use amg_s_krm_solver, amg_protect_name => amg_s_krm_solver_apply
  implicit none
  type(psb_desc_type), intent(in)      :: desc_data
  class(amg_s_krm_solver_type), intent(inout) :: sv
  real(psb_spk_),intent(inout)         :: x(:)
  real(psb_spk_),intent(inout)         :: y(:)
  real(psb_spk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)           :: trans
  real(psb_spk_),target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(out)        :: info
  character, intent(in), optional       :: init
  real(psb_spk_),intent(inout), optional :: initu(:)
  real(psb_spk_), allocatable    :: z(:)
  integer(psb_ipk_) :: np,me,i, err_act, debug_unit, debug_level
  type(psb_ctxt_type) :: ctxt
  character(len=20) :: name='@Z@_krm_solver_apply', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_data%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  call psb_geasb(z,desc_data,info,scratch=.true.)

  call sv%prec%apply(x,z,desc_data,info,trans=trans,work=work)

  call psb_geaxpby(alpha,z,beta,y,desc_data,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine amg_s_krm_solver_apply

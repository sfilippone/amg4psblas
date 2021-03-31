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
! File: amg_cprecbld.f90
!
! Subroutine: amg_cprecbld
! Version:    complex
! Contains:   subroutine init_baseprec_av
!
!  This routine builds the preconditioner according to the requirements made by
!  the user through the subroutines amg_precinit and amg_precset. 
! 
!
! Arguments:
!    a       -  type(psb_cspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_cprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!  
subroutine amg_cprecbld(a,desc_a,prec,info,amold,vmold,imold)

  use psb_base_mod
  use amg_c_prec_mod, amg_protect_name => amg_cprecbld

  Implicit None

  ! Arguments
  type(psb_cspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  class(amg_cprec_type),intent(inout), target         :: prec
  integer(psb_ipk_), intent(out)                               :: info
  class(psb_c_base_sparse_mat), intent(in), optional :: amold
  class(psb_c_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  ! Local Variables
  type(amg_cprec_type) :: t_prec
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  integer(psb_ipk_)   :: err,i,k,err_act, iszv, newsz
  integer(psb_ipk_)   :: ipv(amg_ifpsz_), val
  type(amg_dml_parms) :: prm
  integer(psb_ipk_)   :: debug_level, debug_unit
  character(len=20)   :: name, ch_err

  info=psb_success_
  err=0
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  name = 'amg_cprecbld'
  info = psb_success_
  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  prec%ctxt = ctxt

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '
  !

  if (.not.allocated(prec%precv)) then 
    !! Error: should have called amg_cprecinit
    info=3111
    call psb_errpush(info,name)
    goto 9999
  end if

  !
  ! Check to ensure all procs have the same 
  !   
  newsz = -1
  iszv  = size(prec%precv)
  call psb_bcast(ctxt,iszv)
  if (iszv /= size(prec%precv)) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Inconsistent size of precv')
    goto 9999
  end if

  if (iszv <= 0) then 
    ! Is this really possible? probably not.
    info=psb_err_from_subroutine_
    ch_err='size bpv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  !
  ! Build the  preconditioner
  ! 
  call prec%hierarchy_build(a,desc_a,info)
  
  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Error from hierarchy build')
    goto 9999
  end if
  
  call prec%smoothers_build(a,desc_a,info,amold,vmold,imold)
  
  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Error from smoothers build')
    goto 9999
  end if
  

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_cprecbld

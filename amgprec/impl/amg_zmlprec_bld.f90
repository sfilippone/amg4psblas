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
!  
! File: amg_zmlprec_bld.f90
!
! Subroutine: amg_zmlprec_bld
! Version:    complex
!
!  This routine builds the preconditioner according to the requirements made by
!  the user trough the subroutines amg_precinit and amg_precset.
!  
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level,
!  (for more details see the description of amg_Tonelev_type in amg_prec_type.f90).
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level. No transfer operators are associated to level 1.
!
!  This routine simply calls amg_z_hierarchy_bld and amg_z_smoothers_bld; they
!  can also be called explicitly from the user. 
! 
!
! Arguments:
!    a       -  type(psb_zspmat_type).
!               The sparse matrix structure containing the local part of the
!               matrix to be preconditioned.
!    desc_a  -  type(psb_desc_type), input.
!               The communication descriptor of a.
!    p       -  type(amg_zprec_type), input/output.
!               The preconditioner data structure containing the local part
!               of the preconditioner to be built.
!    info    -  integer, output.
!               Error code.              
!
!    amold   -  class(psb_z_base_sparse_mat), input, optional
!               Mold for the inner format of matrices contained in the
!               preconditioner
!
!
!    vmold   -  class(psb_z_base_vect_type), input, optional
!               Mold for the inner format of vectors contained in the
!               preconditioner
!
!
!  
subroutine amg_zmlprec_bld(a,desc_a,p,info,amold,vmold,imold)

  use psb_base_mod
  use amg_z_inner_mod, amg_protect_name => amg_zmlprec_bld
  use amg_z_prec_mod

  Implicit None

  ! Arguments
  type(psb_zspmat_type),intent(in), target           :: a
  type(psb_desc_type), intent(inout), target           :: desc_a
  type(amg_zprec_type),intent(inout),target          :: p
  integer(psb_ipk_), intent(out)                       :: info
  class(psb_z_base_sparse_mat), intent(in), optional :: amold
  class(psb_z_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  ! Local Variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  integer(psb_ipk_)   :: err,i,k, err_act, iszv, newsz, casize, nplevs, mxplevs
  real(psb_dpk_)       :: mnaggratio
  integer(psb_ipk_)   :: ipv(amg_ifpsz_), val
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

  name = 'amg_zmlprec_bld'
  info = psb_success_
  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Entering '

  call p%hierarchy_build(a,desc_a,info)
  
  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Error from hierarchy build')
    goto 9999
  end if
  
  iszv = p%get_nlevs()

  call p%smoothers_build(a,desc_a,info,amold,vmold,imold)

  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    call psb_errpush(info,name,a_err='Error from smoothers build')
    goto 9999
  end if


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'Exiting with',iszv,' levels'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_zmlprec_bld

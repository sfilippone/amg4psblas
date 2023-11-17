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
subroutine amg_d_poly_smoother_bld(a,desc_a,sm,info,amold,vmold,imold)

  use psb_base_mod
  use amg_d_diag_solver
  use amg_d_l1_diag_solver
  use amg_d_poly_coeff_mod
  use amg_d_poly_smoother, amg_protect_name => amg_d_poly_smoother_bld
  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target          :: a
  Type(psb_desc_type), Intent(inout)                 :: desc_a
  class(amg_d_poly_smoother_type), intent(inout)      :: sm
  integer(psb_ipk_), intent(out)                     :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold
  ! Local variables
  type(psb_dspmat_type) :: tmpa
  integer(psb_ipk_)   :: n_row,n_col, nrow_a, nztota, nzeros
  type(psb_ctxt_type) :: ctxt
  real(psb_dpk_), allocatable :: da(:), dsv(:) 
  integer(psb_ipk_)   :: np, me, i, err_act, debug_unit, debug_level
  character(len=20)   :: name='d_poly_smoother_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  n_row  = desc_a%get_local_rows()
  n_col  = desc_a%get_local_cols()
  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()

  if ((1<=sm%pdegree).and.(sm%pdegree<=30)) then
    call psb_realloc(sm%pdegree,sm%poly_beta,info)
    sm%poly_beta(1:sm%pdegree) = amg_d_beta_mat(1:sm%pdegree,sm%pdegree)
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='invalid sm%degree')
    goto 9999
  end if
  sm%pa => a
  if (.not.allocated(sm%sv)) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,&
         & a_err='unallocated sm%sv')
    goto 9999
  end if
  call sm%sv%build(a,desc_a,info,amold=amold,vmold=vmold)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,&
         & a_err='sv%build')
    goto 9999
  end if
  
  if (.true.) then 
    select type(ssv => sm%sv)
    class is(amg_d_l1_diag_solver_type)
      da  = a%arwsum(info)
      dsv = ssv%dv%get_vect()
      sm%rho_ba = maxval(da(1:n_row)*dsv(1:n_row))
    class default
      write(0,*) 'PolySmoother BUILD: only L1-Jacobi/L1-DIAG for now ',ssv%get_fmt()
      sm%rho_ba = done          
    end select
  else
    block
      type(psb_d_vect_type) :: tq, tz,wv(2)
      real(psb_dpk_)        :: qnrm, lambda
      real(psb_dpk_),allocatable :: work(:)
      integer(psb_ipk_)     :: i, n_cols
      n_cols = desc_a%get_local_cols()
      allocate(work(4*n_cols))
      call psb_geasb(tz,desc_a,info,mold=vmold,scratch=.true.)
      call psb_geasb(wv(1),desc_a,info,mold=vmold,scratch=.true.)
      call psb_geasb(wv(2),desc_a,info,mold=vmold,scratch=.true.)
      call psb_geall(tq,desc_a,info)
      call tq%set(done)
      call psb_geasb(tq,desc_a,info,mold=vmold) 
      call psb_spmm(done,a,tq,dzero,tz,desc_a,info) ! z_1 = A q_0
      do i=1,10
        call sm%sv%apply_v(done,tz,dzero,tq,desc_a,'NoTrans',work,wv,info) ! q_k = M^{-1} q_k
        qnrm = psb_genrmi(tq,desc_a,info)               ! qnrm = |q_k|_inf
        call tq%scal((done/qnrm))                       ! q_k = q_k/qnrm
        call psb_spmm(done,a,tq,dzero,tz,desc_a,info) ! z_{k=1} = A q_k
        lambda = psb_gedot(tq,tz,desc_a,info)      ! lambda = q_k^T z_{k+1} = q_k^T A q_k
        write(0,*) 'BLD: lambda estimate ',i,lambda
      end do
      sm%rho_ba = lambda
    end block
  end if
  
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine amg_d_poly_smoother_bld

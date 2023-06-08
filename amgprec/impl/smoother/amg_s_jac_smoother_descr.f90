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
subroutine amg_s_jac_smoother_descr(sm,info,iout,coarse,prefix)

  use psb_base_mod
  use amg_s_diag_solver
  use amg_s_jac_smoother, amg_protect_name => amg_s_jac_smoother_descr
  use amg_s_diag_solver
  use amg_s_gs_solver

  Implicit None

  ! Arguments
  class(amg_s_jac_smoother_type), intent(in) :: sm
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), intent(in), optional      :: iout
  logical, intent(in), optional              :: coarse
  character(len=*), intent(in), optional  :: prefix

  ! Local variables
  integer(psb_ipk_)      :: err_act
  character(len=20), parameter :: name='amg_s_jac_smoother_descr'
  integer(psb_ipk_)      :: iout_
  logical      :: coarse_
  character(1024)    :: prefix_

  call psb_erractionsave(err_act)
  info = psb_success_
  if (present(coarse)) then 
    coarse_ = coarse
  else
    coarse_ = .false.
  end if
  if (present(iout)) then 
    iout_ = iout 
  else
    iout_ = psb_out_unit
  endif
  if (present(prefix)) then
    prefix_ = prefix
  else
    prefix_ = ""
  end if

  if (.not.coarse_) then
    if (allocated(sm%sv)) then
      select type(smv=>sm%sv)
      class is (amg_s_diag_solver_type)
        write(iout_,*) trim(prefix_), '  Point Jacobi  '
        write(iout_,*) trim(prefix_), '       Local diagonal:'
        call smv%descr(info,iout_,coarse=coarse,prefix=prefix)        
      class is (amg_s_bwgs_solver_type)
        write(iout_,*) trim(prefix_), '  Hybrid Backward Gauss-Seidel  '
      class is (amg_s_gs_solver_type)
        write(iout_,*) trim(prefix_), '  Hybrid Forward Gauss-Seidel  '
      class default
        write(iout_,*) trim(prefix_), '  Block Jacobi  '
        write(iout_,*) trim(prefix_), '       Local solver details:'
        call smv%descr(info,iout_,coarse=coarse,prefix=prefix)
      end select
      
    else
      write(iout_,*) trim(prefix_), '  Block Jacobi  '
    end if
  else
    if (allocated(sm%sv)) then 
      call sm%sv%descr(info,iout_,coarse=coarse,prefix=prefix)
    end if
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine amg_s_jac_smoother_descr

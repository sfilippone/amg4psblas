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
! File: amg_dfile_prec_memory_use.f90
!
!
! Subroutine: amg_file_prec_memory_use
! Version: real
!
!  This routine prints a memory_useiption of the preconditioner to the standard 
!  output or to a file. It must be called after the preconditioner has been
!  built by amg_precbld.
!
! Arguments:
!  p       -  type(amg_Tprec_type), input.
!             The preconditioner data structure to be printed out.
!  info    -  integer, output.
!             error code.
!  iout    -  integer, input, optional.
!             The id of the file where the preconditioner description
!             will be printed. If iout is not present, then the standard
!             output is condidered.
!  root    -  integer, input, optional.
!             The id of the process printing the message; -1 acts as a wildcard.
!             Default is psb_root_
!
!
!
! verbosity:
!        <0: suppress all messages
!         0: normal
!        >1: increased details 
!
subroutine amg_sfile_prec_memory_use(prec,info,iout,root, verbosity,prefix)
  use psb_base_mod
  use amg_s_prec_mod, amg_protect_name => amg_sfile_prec_memory_use
  use amg_s_inner_mod
  use amg_s_gs_solver

  implicit none 
  ! Arguments
  class(amg_sprec_type), intent(in)     :: prec
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), intent(in), optional :: iout
  integer(psb_ipk_), intent(in), optional :: root
  integer(psb_ipk_), intent(in), optional :: verbosity
  character(len=*), intent(in), optional  :: prefix


  ! Local variables
  integer(psb_ipk_)   :: ilev, nlev, ilmin, nswps
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  logical             :: is_symgs
  character(len=20), parameter :: name='amg_file_prec_memory_use'
  integer(psb_ipk_)  :: iout_, root_, verbosity_
  character(1024)    :: prefix_

  info = psb_success_
  if (present(iout)) then 
    iout_ = iout
  else
    iout_ = psb_out_unit
  end if
  if (iout_ < 0) iout_ = psb_out_unit
  if (present(verbosity)) then
    verbosity_ = verbosity
  else
    verbosity_ = 0
  end if
  if (verbosity_ < 0) goto 9998
  if (present(prefix)) then
    prefix_ = prefix
  else
    prefix_ = ""
  end if

  ctxt = prec%ctxt

  if (allocated(prec%precv)) then

    call psb_info(ctxt,me,np)
    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    end if
    if (root_ == -1) root_ = me

    if (verbosity_ >=0) then 
      !
      ! The preconditioner description is printed by processor psb_root_.
      ! This agrees with the fact that all the parameters defining the
      ! preconditioner have the same values on all the procs (this is
      ! ensured by amg_precbld).
      !
      if (me == root_) then

        write(iout_,*) 
        write(iout_,'(a,1x,a)') trim(prefix_),'Preconditioner memory usage'
        nlev = size(prec%precv)
        do ilev=1,nlev
          call prec%precv(ilev)%memory_use(ilev,nlev,ilmin,info, &
               & iout=iout_,verbosity=verbosity,prefix=prefix)
        end do
      end if
    end if
  else
    write(iout_,*) trim(name), &
         & ': Error: no base preconditioner available, something is wrong!'
    info = -2
    return
  endif
9998 continue
end subroutine amg_sfile_prec_memory_use

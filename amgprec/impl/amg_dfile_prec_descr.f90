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
! File: amg_dfile_prec_descr.f90
!
!
! Subroutine: amg_file_prec_descr
! Version: real
!
!  This routine prints a description of the preconditioner to the standard 
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
subroutine amg_dfile_prec_descr(prec,iout,root)
  use psb_base_mod
  use amg_d_prec_mod, amg_protect_name => amg_dfile_prec_descr
  use amg_d_inner_mod
  use amg_d_gs_solver
  
  implicit none 
  ! Arguments
  class(amg_dprec_type), intent(in)      :: prec
  integer(psb_ipk_), intent(in), optional :: iout
  integer(psb_ipk_), intent(in), optional :: root

  ! Local variables
  integer(psb_ipk_)   :: ilev, nlev, ilmin, info, nswps
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  logical             :: is_symgs
  character(len=20), parameter :: name='amg_file_prec_descr'
  integer(psb_ipk_)  :: iout_
  integer(psb_ipk_)  :: root_

  info = psb_success_
  if (present(iout)) then 
    iout_ = iout
  else
    iout_ = psb_out_unit
  end if
  if (iout_ < 0) iout_ = psb_out_unit

  ctxt = prec%ctxt

  if (allocated(prec%precv)) then

    call psb_info(ctxt,me,np)
    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    end if
    if (root_ == -1) root_ = me

    !
    ! The preconditioner description is printed by processor psb_root_.
    ! This agrees with the fact that all the parameters defining the
    ! preconditioner have the same values on all the procs (this is
    ! ensured by amg_precbld).
    !
    if (me == root_) then
      nlev = size(prec%precv)
      do ilev = 1, nlev 
        if (.not.allocated(prec%precv(ilev)%sm)) then 
          info = 3111
          write(iout_,*) ' ',name,&
               & ': error: inconsistent MLPREC part, should call amg_PRECINIT'
          return
        endif
      end do

      write(iout_,*) 
      write(iout_,'(a)') 'Preconditioner description'

      if (nlev == 1) then
        !
        ! Here we have a gigantic kludge just to handle Symmetrized Gauss-Seidel.
        ! Will need rethinking...
        !
        if (allocated(prec%precv(1)%sm2a)) then
          is_symgs = .false.
          select type(sv2 => prec%precv(1)%sm2a%sv)
          class is (amg_d_bwgs_solver_type)
            select type(sv1 => prec%precv(1)%sm%sv)
            class is (amg_d_gs_solver_type)
              is_symgs = .true.
            end select
          end select
          if (is_symgs) then
            write(iout_,*) ' Forward-Backward (symmetrized) Hybrid Gauss-Seidel'
          else
            write(iout_,*) 'Pre Smoother: '
            call prec%precv(1)%sm%descr(info,iout=iout_)
            write(iout_,*) 'Post smoother:'
            call prec%precv(1)%sm2a%descr(info,iout=iout_)
          end if
          nswps = max(prec%precv(1)%parms%sweeps_pre,prec%precv(1)%parms%sweeps_post)

        else
          call prec%precv(1)%sm%descr(info,iout=iout_)
          nswps = prec%precv(1)%parms%sweeps_pre
        end if
        if (nswps > 1)  write(iout_,*) '  Number of  sweeps : ',nswps
        write(iout_,*) 

      else if (nlev > 1) then
        !
        ! Print description of base preconditioner
        !
        write(iout_,*) 'Multilevel Preconditioner'
        write(iout_,*) 'Outer sweeps:',prec%outer_sweeps
        write(iout_,*) 
        if (allocated(prec%precv(1)%sm2a)) then
          write(iout_,*) 'Pre Smoother: '
          call prec%precv(1)%sm%descr(info,iout=iout_)
          write(iout_,*) 'Post smoother:'
          call prec%precv(1)%sm2a%descr(info,iout=iout_)
        else
          write(iout_,*) 'Smoother: '
          call prec%precv(1)%sm%descr(info,iout=iout_)
        end if
        !
        ! Print multilevel details
        !
        write(iout_,*) 
        write(iout_,*) 'Multilevel hierarchy: '
        write(iout_,*) ' Number of levels   : ',nlev
        write(iout_,*) ' Operator complexity: ',prec%get_complexity()
        write(iout_,*) ' Average coarsening : ',prec%get_avg_cr()
        ilmin = 2
        if (nlev == 2) ilmin=1
        do ilev=ilmin,nlev
          call prec%precv(ilev)%descr(ilev,nlev,ilmin,info,iout=iout_)
        end do
        write(iout_,*) 

      else
        write(iout_,*) trim(name), &
             & ': invalid preconditioner array size ?',nlev
        info = -2
        return

      end if
    end if

  else
    write(iout_,*) trim(name), &
         & ': Error: no base preconditioner available, something is wrong!'
    info = -2
    return
  endif

end subroutine amg_dfile_prec_descr

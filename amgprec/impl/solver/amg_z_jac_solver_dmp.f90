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
subroutine amg_z_jac_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)

  use psb_base_mod
  use amg_z_jac_solver, amg_protect_name => amg_z_jac_solver_dmp
  implicit none 
  class(amg_z_jac_solver_type), intent(in) :: sv
  type(psb_desc_type), intent(in)             :: desc
  integer(psb_ipk_), intent(in)              :: level
  integer(psb_ipk_), intent(out)             :: info
  character(len=*), intent(in), optional     :: prefix, head
  logical, optional, intent(in)              :: solver, global_num
  integer(psb_ipk_)   :: i, j, il1, iln, lname, lev
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np
  character(len=80)   :: prefix_
  character(len=120)  :: fname ! len should be at least 20 more than
  logical :: solver_, global_num_
  integer(psb_lpk_), allocatable :: iv(:)
  !  len of prefix_ 

  info = 0

  ctxt = desc%get_context()
  call psb_info(ctxt,iam,np)

  if (present(solver)) then 
    solver_ = solver
  else
    solver_ = .false. 
  end if
  if (present(global_num)) then 
    global_num_ = global_num
  else
    global_num_ = .false. 
  end if


  if (solver_) then 
    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_slv_z"
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5
    
    if (global_num_) then
      iv = desc%get_global_indices(owned=.false.)      

      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_nd.mtx'
      if (sv%a%is_asb()) &
           & call sv%a%print(fname,head=head,iv=iv)
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_diag.mtx'
      if (allocated(sv%dv)) &
           & call psb_geprt(fname,sv%dv%v%v,head=head)
      
    else

      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_nd.mtx'
      if (sv%a%is_asb()) &
           & call sv%a%print(fname,head=head)
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_diag.mtx'
      if (allocated(sv%dv)) &
           & call psb_geprt(fname,sv%dv%v%v,head=head)
    end if
  end if

end subroutine amg_z_jac_solver_dmp

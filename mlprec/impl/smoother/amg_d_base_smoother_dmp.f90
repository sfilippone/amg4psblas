!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
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
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine amg_d_base_smoother_dmp(sm,desc,level,info,prefix,head,smoother,solver,global_num)
  
  use psb_base_mod
  use amg_d_base_smoother_mod, amg_protect_name =>  amg_d_base_smoother_dmp
  implicit none 
  class(amg_d_base_smoother_type), intent(in) :: sm
  type(psb_desc_type), intent(in)               :: desc
  integer(psb_ipk_), intent(in)              :: level
  integer(psb_ipk_), intent(out)             :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: smoother, solver, global_num
  integer(psb_ipk_)  :: i, j, il1, iln, lname, lev
  integer(psb_ipk_)  :: ictxt,iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than
  logical :: smoother_, global_num_
  !  len of prefix_ 

  info = 0

  if (present(prefix)) then 
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_smth_d"
  end if
  ictxt = desc%get_context()
  call psb_info(ictxt,iam,np)

  if (present(global_num)) then 
    global_num_ = global_num
  else
    global_num_ = .false. 
  end if
  if (present(smoother)) then 
    smoother_ = smoother
  else
    smoother_ = .false. 
  end if
  lname = len_trim(prefix_)
  fname = trim(prefix_)
  write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
  lname = lname + 5

  ! At base level do nothing for the smoother
  if (allocated(sm%sv)) &
       & call sm%sv%dump(desc,level,info,solver=solver,prefix=prefix,global_num=global_num)

end subroutine amg_d_base_smoother_dmp

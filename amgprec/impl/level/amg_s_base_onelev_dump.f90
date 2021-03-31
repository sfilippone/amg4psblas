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
subroutine amg_s_base_onelev_dump(lv,level,info,prefix,head,ac,rp,&
     & smoother,solver,tprol,global_num)

  use psb_base_mod
  use amg_s_onelev_mod, amg_protect_name => amg_s_base_onelev_dump
  implicit none 
  class(amg_s_onelev_type), intent(in) :: lv
  integer(psb_ipk_), intent(in)          :: level
  integer(psb_ipk_), intent(out)         :: info
  character(len=*), intent(in), optional :: prefix, head
  logical, optional, intent(in)    :: ac, rp, smoother, solver, tprol, global_num
  ! Local variables
  integer(psb_ipk_)   :: i, j, il1, iln, lname, lev, ni
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np
  character(len=80)   :: prefix_, frmt
  character(len=1024) :: fname 
  logical :: ac_, rp_, tprol_, global_num_
  integer(psb_lpk_), allocatable :: ivr(:), ivc(:)

  info = 0

  if (present(prefix)) then 
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_lev_s"
  end if

  if (associated(lv%base_desc)) then 
    ctxt = lv%base_desc%get_context()
    call psb_info(ctxt,iam,np)
  else 
    iam     = -1
    np      = -1
  end if
  if (present(ac)) then 
    ac_ = ac
  else
    ac_ = .false. 
  end if
  if (present(rp)) then 
    rp_ = rp
  else
    rp_ = .false. 
  end if
  if (present(tprol)) then 
    tprol_ = tprol
  else
    tprol_ = .false. 
  end if
  if (present(global_num)) then 
    global_num_ = global_num
  else
    global_num_ = .false. 
  end if
  lname = len_trim(prefix_)
  fname = trim(prefix_)

  if (np > 0) then 
    ni  = floor(log10(1.0*np)) + 1
    write(frmt,'(a,i3.3,a,i3.3,a)') '(a,i',ni,'.',ni,')'
    write(fname(lname+1:lname+ni+2),frmt) '_p',iam
    lname = lname + ni + 2
  end if

  if (global_num_) then
    if (level >= 2) then 
      if (ac_) then
        ivr = lv%desc_ac%get_global_indices(owned=.false.)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
        call lv%ac%print(fname,head=head,iv=ivr)
      end if
      if (rp_) then 
        ivr = lv%linmap%p_desc_U%get_global_indices(owned=.false.)
        ivc = lv%linmap%p_desc_V%get_global_indices(owned=.false.)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_r.mtx'
        call lv%linmap%mat_U2V%print(fname,head=head,ivr=ivc,ivc=ivr)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_p.mtx'
        call lv%linmap%mat_V2U%print(fname,head=head,ivr=ivr,ivc=ivc)
      end if
      if (tprol_) then
        ! Tentative prolongator is stored with column indices already
        ! in global numbering, so only IVR is needed. 
        ivr = lv%linmap%p_desc_U%get_global_indices(owned=.false.)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_tprol.mtx'
        !
        !  This is not implemented yet.
        !call lv%tprol%print(fname,head=head,ivr=ivr)
      end if
    end if
  else
    if (level >= 2) then 
      if (ac_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_ac.mtx'
        call lv%ac%print(fname,head=head)
      end if
      if (rp_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_r.mtx'
        call lv%linmap%mat_U2V%print(fname,head=head)
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_p.mtx'
        call lv%linmap%mat_V2U%print(fname,head=head)
      end if
      if (tprol_) then 
        write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_tprol.mtx'
        !
        !  This is not implemented yet.
        !call lv%tprol%print(fname,head=head)
      end if
    end if
  end if
  
  if (level >= 2) then 
    if (allocated(lv%sm)) then
      call lv%sm%dump(lv%desc_ac,level,info,smoother=smoother, &
           & solver=solver,prefix=trim(prefix_)//"_sm",global_num=global_num)
    end if
    if (allocated(lv%sm2a)) then
      call lv%sm2a%dump(lv%desc_ac,level,info,smoother=smoother, &
           & solver=solver,prefix=trim(prefix_)//"_sm2a",global_num=global_num)
    end if
  else
    if (allocated(lv%sm)) then
      call lv%sm%dump(lv%base_desc,level,info,smoother=smoother, &
           & solver=solver,prefix=trim(prefix_)//"_sm",global_num=global_num)
    end if
    if (allocated(lv%sm2a)) then
      call lv%sm2a%dump(lv%base_desc,level,info,smoother=smoother, &
           & solver=solver,prefix=trim(prefix_)//"_sm2a",global_num=global_num)
    end if
  end if
  
end subroutine amg_s_base_onelev_dump

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
subroutine amg_s_base_onelev_setsv(lev,val,info,pos)

  use psb_base_mod
  use amg_s_prec_mod, amg_protect_name => amg_s_base_onelev_setsv

  implicit none

  ! Arguments
  class(amg_s_onelev_type), target, intent(inout) :: lev
  class(amg_s_base_solver_type), intent(in)       :: val
  integer(psb_ipk_), intent(out)                  :: info
  character(len=*), optional, intent(in)          :: pos
  
  ! Local variables
  integer(psb_ipk_)                :: ipos_
  character(len=*), parameter      :: name='amg_base_onelev_setsv'

  info = psb_success_

  if (present(pos)) then
    select case(psb_toupper(trim(pos)))
    case('PRE')
      ipos_ = amg_smooth_pre_
    case('POST')
      ipos_ = amg_smooth_post_
    case default
      ipos_ = amg_smooth_both_
    end select
  else
    ipos_ = amg_smooth_both_
  end if
  
  if ((ipos_ == amg_smooth_pre_).or.(ipos_ == amg_smooth_both_)) then 
    if (allocated(lev%sm)) then 
      if (allocated(lev%sm%sv)) then
        if (.not.same_type_as(lev%sm%sv,val))  then
          call lev%sm%sv%free(info)
          if (info == 0) deallocate(lev%sm%sv,stat=info)
          if (info /= 0) then
            info = 3111
            return
          end if
        end if
      end if
      
      if (.not.allocated(lev%sm%sv)) then 
#ifdef HAVE_MOLD 
        allocate(lev%sm%sv,mold=val,stat=info) 
#else
        allocate(lev%sm%sv,source=val,stat=info) 
#endif
        if (info /= 0) then
          info = 3111
          return
        end if
      end if
      call lev%sm%sv%default()
    else
      info = 3111
      write(psb_err_unit,*) name,&
           &': Error: uninitialized preconditioner component,',&
           &' should call amg_PRECINIT/amg_PRECSET' 
      return 
      
    end if
  end if

  !
  ! If POS was not specified and therefore we have amg_smooth_both_
  ! we need to update sm2a *only* if it was already allocated,
  ! otherwise it is not needed (since we have just fixed %sm in the
  ! pre section). 
  !

  if ((ipos_ == amg_smooth_post_).or. &
       ((ipos_ == amg_smooth_both_).and.(allocated(lev%sm2a)))) then 


    if (allocated(lev%sm2a)) then 
      if (allocated(lev%sm2a%sv)) then
        if (.not.same_type_as(lev%sm2a%sv,val))  then
          call lev%sm2a%sv%free(info)
          if (info == 0) deallocate(lev%sm2a%sv,stat=info)
          if (info /= 0) then
            info = 3111
            return
          end if
        end if
      end if
      if (.not.allocated(lev%sm2a%sv)) then 
#ifdef HAVE_MOLD 
        allocate(lev%sm2a%sv,mold=val,stat=info) 
#else
        allocate(lev%sm2a%sv,source=val,stat=info) 
#endif
        if (info /= 0) then
          info = 3111
          return
        end if
      end if
      call lev%sm2a%sv%default()
      
    else
      info = 3111
      write(psb_err_unit,*) name,&
           &': Error: uninitialized preconditioner component,',&
           &' should call amg_PRECINIT/amg_PRECSET' 
      return 
      
    end if
    
  end if
  
end subroutine amg_s_base_onelev_setsv


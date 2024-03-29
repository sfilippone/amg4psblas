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
! File: amg_cprecset.f90
!
subroutine amg_cprecsetsm(p,val,info,ilev,ilmax,pos)

  use psb_base_mod
  use amg_c_prec_mod, amg_protect_name => amg_cprecsetsm

  implicit none

  ! Arguments
  class(amg_cprec_type), target, intent(inout):: p
  class(amg_c_base_smoother_type), intent(in) :: val
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
  character(len=*), optional, intent(in)      :: pos

  ! Local variables
  integer(psb_ipk_)                      :: ilev_, nlev_, ilmin_, ilmax_
  character(len=*), parameter            :: name='amg_precsetsm'

  info = psb_success_

  if (.not.allocated(p%precv)) then 
    info = 3111
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call amg_PRECINIT'
    return 
  endif
  nlev_ = size(p%precv)

  if (present(ilev)) then 
    ilev_ = ilev
    ilmin_ = ilev
    if (present(ilmax)) then
      ilmax_ = ilmax
    else
      ilmax_ = ilev_
    end if
  else
    ilev_ = 1 
    ilmin_ = 1
    ilmax_ = nlev_
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if ((ilmax_<1).or.(ilmax_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILMAX/NLEV combination',ilmax_, nlev_
    return
  endif  

  do ilev_ = ilmin_, ilmax_ 
    call p%precv(ilev_)%set(val,info,pos=pos)
    if (info /= 0) return 
  end do

end subroutine amg_cprecsetsm

subroutine amg_cprecsetsv(p,val,info,ilev,ilmax,pos)

  use psb_base_mod
  use amg_c_prec_mod, amg_protect_name => amg_cprecsetsv

  implicit none

  ! Arguments
  class(amg_cprec_type), intent(inout)       :: p
  class(amg_c_base_solver_type), intent(in) :: val
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), optional, intent(in)   :: ilev,ilmax
  character(len=*), optional, intent(in)      :: pos

  ! Local variables
  integer(psb_ipk_)                       :: ilev_, nlev_, ilmin_, ilmax_
  character(len=*), parameter            :: name='amg_precsetsv'

  info = psb_success_

  if (.not.allocated(p%precv)) then 
    info = 3111
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call amg_PRECINIT'
    return 
  endif
  nlev_ = size(p%precv)

  if (present(ilev)) then 
    ilev_ = ilev
    ilmin_ = ilev
    if (present(ilmax)) then
      ilmax_ = ilmax
    else
      ilmax_ = ilev_
    end if
  else
    ilev_ = 1 
    ilmin_ = 1
    ilmax_ = nlev_
  end if

  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if ((ilmax_<1).or.(ilmax_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILMAX/NLEV combination',ilmax_, nlev_
    return
  endif  

  do ilev_ = ilmin_, ilmax_ 
    call p%precv(ilev_)%set(val,info,pos=pos)
    if (info /= 0) return 
  end do

end subroutine amg_cprecsetsv

subroutine amg_cprecsetag(p,val,info,ilev,ilmax,pos)

  use psb_base_mod
  use amg_c_prec_mod, amg_protect_name => amg_cprecsetag

  implicit none

  ! Arguments
  class(amg_cprec_type), intent(inout)       :: p
  class(amg_c_base_aggregator_type), intent(in) :: val
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), optional, intent(in)   :: ilev, ilmax
  character(len=*), optional, intent(in)      :: pos

  ! Local variables
  integer(psb_ipk_)                       :: ilev_, nlev_, ilmin_, ilmax_
  character(len=*), parameter            :: name='amg_precsetag'

  info = psb_success_

  if (.not.allocated(p%precv)) then 
    info = 3111
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call amg_PRECINIT'
    return 
  endif
  nlev_ = size(p%precv)

  if (present(ilev)) then 
    ilev_ = ilev
    ilmin_ = ilev
    if (present(ilmax)) then
      ilmax_ = ilmax
    else
      ilmax_ = ilev_
    end if
  else
    ilev_ = 1 
    ilmin_ = 1
    ilmax_ = nlev_
  end if


  if ((ilev_<1).or.(ilev_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         & ': Error: invalid ILEV/NLEV combination',ilev_, nlev_
    return
  endif
  if ((ilmax_<1).or.(ilmax_ > nlev_)) then 
    info = -1
    write(psb_err_unit,*) name,&
         &': Error: invalid ILMAX/NLEV combination',ilmax_, nlev_
    return
  endif  

  do ilev_ = ilmin_, ilmax_ 
    call p%precv(ilev_)%set(val,info,pos=pos)
    if (info /= 0) return 
  end do

end subroutine amg_cprecsetag


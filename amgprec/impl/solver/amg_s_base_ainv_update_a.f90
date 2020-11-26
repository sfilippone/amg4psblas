!
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
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
subroutine amg_s_base_ainv_update_a(sv,x,desc_data,info)
  use amg_s_base_ainv_mod, amg_protect_name => amg_s_base_ainv_update_a
  use psb_base_mod

  implicit none

  type(psb_desc_type), intent(in)            :: desc_data
  class(amg_s_base_ainv_solver_type), intent(inout) :: sv
  real(psb_spk_),intent(in)                  :: x(:)
  integer(psb_ipk_), intent(out)             :: info

  ! Local variables
  real(psb_spk_), allocatable   :: dd(:), ee(:)
  integer(psb_ipk_) :: nrows, ncols, i, j, k, nzr, nzc, ir, ic

  type(psb_s_csc_sparse_mat) :: ac
  type(psb_s_csr_sparse_mat) :: ar

  dd = sv%dv%get_vect()
  if (size(x)<size(dd)) then
    write(0,*) 'Wrong size in update',size(x),size(dd)
    info = -1
    return
  end if

  ee = x(1:size(dd))
  call sv%z%cp_to(ac)
  call sv%w%cp_to(ar)
  ! Compute diag(W ee Z)
  ! First  ee Z
  call ac%scal(ee,info,side='L')
  ! Then  diag(W (ee Z))
  do i=1,ar%get_nrows()
    ir  = ar%irp(i);  nzr = ar%irp(i+1)-ar%irp(i)
    ic  = ac%icp(i);  nzc = ac%icp(i+1)-ac%icp(i)
    ee(i) = psb_spdot_srtd(&
         & nzr,ar%ja(ir:ir+nzr-1),ar%val(ir:ir+nzr-1),&
         & nzc,ac%ia(ic:ic+nzc-1),ac%val(ic:ic+nzc-1))
  end do
  ! Now invert diagonal
  dd = done/dd
  dd = dd + ee
  dd = done/dd
  sv%d = dd
  call sv%dv%bld(dd,mold=sv%dv%v)

end subroutine amg_s_base_ainv_update_a

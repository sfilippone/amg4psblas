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
!
module amg_base_ainv_mod

  use amg_base_prec_type

  integer, parameter   :: amg_inv_fillin_ = amg_ifpsz_ + 1
  integer, parameter   :: amg_ainv_alg_   = amg_inv_fillin_ + 1
  integer, parameter   :: amg_inv_thresh_ = amg_max_sub_solve_ + 1
#if 0
  integer, parameter   :: amg_ainv_orth1_ = amg_inv_thresh_ + 1
  integer, parameter   :: amg_ainv_orth2_ = amg_ainv_orth1_ + 1
  integer, parameter   :: amg_ainv_orth3_ = amg_ainv_orth2_ + 1
  integer, parameter   :: amg_ainv_orth4_ = amg_ainv_orth3_ + 1
  integer, parameter   :: amg_ainv_llk_   = amg_ainv_orth4_ + 1
#else
  integer, parameter   :: amg_ainv_llk_   = amg_inv_thresh_ + 1
#endif
  integer, parameter   :: amg_ainv_s_llk_ = amg_ainv_llk_ + 1
  integer, parameter   :: amg_ainv_s_ft_llk_ = amg_ainv_s_llk_ + 1
  integer, parameter   :: amg_ainv_llk_noth_ = amg_ainv_s_ft_llk_  + 1
  integer, parameter   :: amg_ainv_mlk_      = amg_ainv_llk_noth_  + 1
  integer, parameter   :: amg_ainv_lmx_      = amg_ainv_mlk_
#if defined(HAVE_TUMA_SAINV)
  integer, parameter   :: amg_ainv_s_tuma_    = amg_ainv_lmx_  + 1
  integer, parameter   :: amg_ainv_l_tuma_    = amg_ainv_s_tuma_  + 1
#endif


end module amg_base_ainv_mod

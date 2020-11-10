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
subroutine amg_s_base_onelev_map_prol_v(lv,alpha,vect_v,beta,vect_u,info,work,vtx,vty)
  use psb_base_mod
  use amg_s_onelev_mod, amg_protect_name =>  amg_s_base_onelev_map_prol_v
  
  implicit none
  class(amg_s_onelev_type), target, intent(inout) :: lv
  real(psb_spk_), intent(in)           :: alpha, beta
  type(psb_s_vect_type), intent(inout) :: vect_u, vect_v
  integer(psb_ipk_), intent(out)       :: info
  real(psb_spk_), optional          :: work(:)
  type(psb_s_vect_type), optional, target, intent(inout)  :: vtx,vty
  
  if (lv%remap_data%ac_pre_remap%is_asb()) then
    !
    ! Remap has happened, deal with it
    !
    write(0,*) 'Remap handling not implemented yet '
  else
    ! Default transfer
    call lv%linmap%map_V2U(alpha,vect_v,beta,vect_u,info,&
         & work=work,vtx=vtx,vty=vty)
  end if
  
end subroutine amg_s_base_onelev_map_prol_v

subroutine amg_s_base_onelev_map_prol_a(lv,alpha,v,beta,u,info,work)
  use psb_base_mod
  use amg_s_onelev_mod, amg_protect_name =>  amg_s_base_onelev_map_prol_a
  implicit none
  class(amg_s_onelev_type), target, intent(inout) :: lv
  real(psb_spk_), intent(in)     :: alpha, beta
  real(psb_spk_), intent(inout)   :: u(:)
  real(psb_spk_), intent(out)     :: v(:)
  integer(psb_ipk_), intent(out) :: info
  real(psb_spk_), optional          :: work(:)

  if (lv%remap_data%ac_pre_remap%is_asb()) then
    !
    ! Remap has happened, deal with it
    !
    write(0,*) 'Remap handling not implemented yet '
  else
    ! Default transfer
    call lv%linmap%map_V2U(alpha,v,beta,u,info,&
         & work=work)
  end if
  
end subroutine amg_s_base_onelev_map_prol_a

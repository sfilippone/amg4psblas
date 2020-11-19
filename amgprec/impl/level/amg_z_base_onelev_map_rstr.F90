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

subroutine amg_z_base_onelev_map_rstr_v(lv,alpha,vect_u,beta,vect_v,info,&
     &  work,vtx,vty)
  use psb_base_mod
  use amg_z_onelev_mod, amg_protect_name =>  amg_z_base_onelev_map_rstr_v
  implicit none
  class(amg_z_onelev_type), target, intent(inout) :: lv
  complex(psb_dpk_), intent(in)           :: alpha, beta
  type(psb_z_vect_type), intent(inout) :: vect_u, vect_v
  integer(psb_ipk_), intent(out)       :: info
  complex(psb_dpk_), optional          :: work(:)
  type(psb_z_vect_type), optional, target, intent(inout)  :: vtx,vty

!!$  write(0,*) 'New map_rstr',lv%remap_data%ac_pre_remap%is_asb()
  if (lv%remap_data%ac_pre_remap%is_asb()) then
    !
    ! Remap has happened, deal with it
    !
!!$    write(0,*) 'Remap handling not implemented yet '
    block
      type(psb_ctxt_type) :: ctxt, nctxt
      integer(psb_ipk_) :: i,j,ip, idest, nsrc, nrl, kp
      integer(psb_ipk_) :: me, np,  rme, rnp
      complex(psb_dpk_), allocatable :: rsnd(:), rrcv(:)
      type(psb_z_vect_type) :: tv
      
      ctxt = lv%remap_data%desc_ac_pre_remap%get_ctxt()
      call psb_info(ctxt,me,np)
      nctxt = lv%desc_ac%get_ctxt()
      call psb_info(nctxt,rme,rnp)
!!$      write(0,*) 'New context ',rme,rnp
      idest = lv%remap_data%idest
      associate(isrc => lv%remap_data%isrc, nrsrc => lv%remap_data%nrsrc)
!!$        write(0,*) 'Should apply maps, then send data from ',me,' to ',idest
!!$        if (rme >= 0) write(0,*) rme, '  Receiving      data from ',isrc(:)
        nsrc = size(isrc)
        nrl  = lv%remap_data%desc_ac_pre_remap%get_local_rows()
        call psb_geall(tv,lv%remap_data%desc_ac_pre_remap,info)
        call psb_geasb(tv,lv%remap_data%desc_ac_pre_remap,info) 
!!$        write(0,*) me,' Size of TV ',tv%get_nrows()
        call lv%linmap%map_U2V(alpha,vect_u,beta,tv,info,&
             & work=work,vtx=vtx,vty=vty)
        rsnd = tv%get_vect()
        call psb_snd(ctxt,rsnd(1:nrl),idest)
        if (rme >=0) then
          allocate(rrcv(sum(nrsrc)))
!!$          write(0,*) me,rme,' Size check ',size(rrcv)!,lv%desc_ac%get_local_rows()
          kp = 0
          do i = 1,size(isrc)
            ip = isrc(i)
            nrl = nrsrc(i)
!!$            write(0,*) me,' Receiving from ',ip,nrl,kp+1,kp+nrl,size(rrcv)
            call psb_rcv(ctxt,rrcv(kp+1:kp+nrl),ip)
            kp = kp + nrl
          end do
          call vect_v%set_vect(rrcv)
        end if
      end associate
!!$      write(0,*) me, ' Restrictor with remap done '
    end block
 
  else
    ! Default transfer
    call lv%linmap%map_U2V(alpha,vect_u,beta,vect_v,info,&
         & work=work,vtx=vtx,vty=vty)
  end if
  
end subroutine amg_z_base_onelev_map_rstr_v

subroutine amg_z_base_onelev_map_rstr_a(lv,alpha,u,beta,v,info,work)
  use psb_base_mod
  use amg_z_onelev_mod, amg_protect_name =>  amg_z_base_onelev_map_rstr_a
  implicit none
  class(amg_z_onelev_type), target, intent(inout) :: lv
  complex(psb_dpk_), intent(in)     :: alpha, beta
  complex(psb_dpk_), intent(inout)   :: u(:)
  complex(psb_dpk_), intent(out)     :: v(:)
  integer(psb_ipk_), intent(out) :: info
  complex(psb_dpk_), optional       :: work(:)

  if (lv%remap_data%ac_pre_remap%is_asb()) then
    !
    ! Remap has happened, deal with it
    !
    write(0,*) 'Remap handling not implemented yet '
  else
    ! Default transfer
    call lv%linmap%map_U2V(alpha,u,beta,v,info,&
         & work=work)
  end if
  
end subroutine amg_z_base_onelev_map_rstr_a

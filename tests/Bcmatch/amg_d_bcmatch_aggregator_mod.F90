!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008, 2010, 2012, 2015, 2017 , 2017 
!  
!                        Salvatore Filippone  Cranfield University
!  		      Ambra Abdullahi Hassan University of Rome Tor Vergata
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
!        Daniela di Serafino    University of Campania "L. Vanvitelli", Caserta, IT
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
!
!
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. This variant is based on the hybrid method
!  presented in 
!
!    S. Gratton, P. Henon, P. Jiranek and X. Vasseur:
!    Reducing complexity of algebraic multigrid by aggregation
!    Numerical Lin. Algebra with Applications, 2016, 23:501-518
!    
  !
  !   sm           -  class(mld_T_base_smoother_type), allocatable
  !                   The current level preconditioner (aka smoother).
  !   parms        -  type(mld_RTml_parms)
  !                   The parameters defining the multilevel strategy.
  !   ac           -  The local part of the current-level matrix, built by
  !                   coarsening the previous-level matrix.
  !   desc_ac      -  type(psb_desc_type).
  !                   The communication descriptor associated to the matrix
  !                   stored in ac.
  !   base_a       -  type(psb_Tspmat_type), pointer.
  !                   Pointer (really a pointer!) to the local part of the current 
  !                   matrix (so we have a unified treatment of residuals).
  !                   We need this to avoid passing explicitly the current matrix
  !                   to the routine which applies the preconditioner.
  !   base_desc    -  type(psb_desc_type), pointer.
  !                   Pointer to the communication descriptor associated to the
  !                   matrix pointed by base_a.
  !   map          -  Stores the maps (restriction and prolongation) between the
  !                   vector spaces associated to the index spaces of the previous
  !                   and current levels.
  !
  !   Methods:  
  !     Most methods follow the encapsulation hierarchy: they take whatever action
  !     is appropriate for the current object, then call the corresponding method for
  !     the contained object.
  !     As an example: the descr() method prints out a description of the
  !     level. It starts by invoking the descr() method of the parms object,
  !     then calls the descr() method of the smoother object. 
  !
  !    descr      -   Prints a description of the object.
  !    default    -   Set default values
  !    dump       -   Dump to file object contents
  !    set        -   Sets various parameters; when a request is unknown
  !                   it is passed to the smoother object for further processing.
  !    check      -   Sanity checks.
  !    sizeof     -   Total memory occupation in bytes
  !    get_nzeros -   Number of nonzeros 
  !
  !

module bcm_csr_type_mod
 use iso_c_binding
 type, bind(c)::  bcm_Vector
   type(c_ptr) :: data
   integer(c_int) :: size
   integer(c_int) :: owns_data
 end type 

 type, bind(c)::  bcm_CSRMatrix
   type(c_ptr) :: i
   type(c_ptr) :: j
   integer(c_int) :: num_rows
   integer(c_int) :: num_cols
   integer(c_int) :: num_nonzeros
   integer(c_int) :: owns_data
   type(c_ptr) :: data
 end type 
end module bcm_csr_type_mod

module mld_d_bcmatch_aggregator_mod
  use mld_d_base_aggregator_mod
  use bcm_csr_type_mod

  type, extends(mld_d_base_aggregator_type) :: mld_d_bcmatch_aggregator_type
    integer(psb_ipk_) :: matching_alg
    integer(psb_ipk_) :: n_sweeps
    real(psb_dpk_), allocatable :: w_tmp(:), w_nxt(:) 
    type(bcm_Vector)  :: w_par
    integer(psb_ipk_) :: max_csize
    integer(psb_ipk_) :: max_nlevels
  contains
    procedure, pass(ag) :: bld_tprol    => mld_d_bcmatch_aggregator_build_tprol
    procedure, pass(ag) :: cseti        => d_bcmatch_aggr_cseti
    procedure, pass(ag) :: default      => d_bcmatch_aggr_set_default
    procedure, pass(ag) :: mat_asb      => mld_d_bcmatch_aggregator_mat_asb
    procedure, pass(ag) :: update_next => d_bcmatch_aggregator_update_next
    procedure, pass(ag) :: bld_wnxt     => d_bcmatch_bld_wnxt
    procedure, pass(ag) :: bld_default_w     => d_bld_default_w
    procedure, pass(ag) :: set_c_default_w     => d_set_default_bcm_w
    procedure, pass(ag) :: descr        => d_bcmatch_aggregator_descr
    procedure, pass(ag) :: clone        => d_bcmatch_aggregator_clone
    procedure, pass(ag) :: free         => d_bcmatch_aggregator_free
!!$    procedure, pass(ag) :: default      => mld_d_base_aggregator_default
    procedure, nopass   :: fmt          => d_bcmatch_aggregator_fmt
  end type mld_d_bcmatch_aggregator_type


  interface
    subroutine  mld_d_bcmatch_aggregator_build_tprol(ag,parms,a,desc_a,ilaggr,nlaggr,op_prol,info)
      import :: mld_d_bcmatch_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)  :: parms 
      type(psb_dspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_dspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine mld_d_bcmatch_aggregator_build_tprol
  end interface

  interface
    subroutine  mld_d_bcmatch_aggregator_mat_asb(ag,parms,a,desc_a,ilaggr,nlaggr,ac,&
         & op_prol,op_restr,info)
      import :: mld_d_bcmatch_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
      type(mld_dml_parms), intent(inout)   :: parms 
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(inout)     :: ilaggr(:), nlaggr(:)
      type(psb_dspmat_type), intent(inout)   :: op_prol
      type(psb_dspmat_type), intent(out)   :: ac,op_restr
      integer(psb_ipk_), intent(out)       :: info
    end subroutine mld_d_bcmatch_aggregator_mat_asb
  end interface  
  

  interface
    subroutine mld_d_bcmatch_map_to_tprol(desc_a,ilaggr,nlaggr,valaggr, op_prol,info)
      import :: mld_d_bcmatch_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), allocatable, intent(inout)  :: ilaggr(:),nlaggr(:)
      real(psb_dpk_), allocatable, intent(inout)  :: valaggr(:)
      type(psb_dspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_d_bcmatch_map_to_tprol
  end interface
  
  interface
    subroutine mld_daggrmat_unsmth_spmm_asb(a,desc_a,ilaggr,nlaggr,parms,ac,op_prol,op_restr,info)
      import :: mld_d_bcmatch_aggregator_type, psb_desc_type, psb_dspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_long_int_k_, mld_dml_parms
      implicit none
      type(psb_dspmat_type), intent(in)        :: a
      type(psb_desc_type), intent(in)            :: desc_a
      integer(psb_ipk_), intent(inout)           :: ilaggr(:), nlaggr(:)
      type(mld_dml_parms), intent(inout)      :: parms 
      type(psb_dspmat_type), intent(inout)     :: op_prol
      type(psb_dspmat_type), intent(out)       :: ac,op_restr
      integer(psb_ipk_), intent(out)             :: info
    end subroutine mld_daggrmat_unsmth_spmm_asb
  end interface
  
  
contains

  subroutine d_bld_default_w(ag,nr)
    use psb_realloc_mod
    implicit none 
    class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_ipk_), intent(in) :: nr
    integer(psb_ipk_) :: info
    call psb_realloc(nr,ag%w_tmp,info)
    if (info /= psb_success_) return
    ag%w_tmp = done
    call ag%set_c_default_w()
  end subroutine d_bld_default_w

  subroutine d_set_default_bcm_w(ag)
    use psb_realloc_mod
    use iso_c_binding
    implicit none 
    class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag

    ag%w_par%size      = psb_size(ag%w_tmp)
    ag%w_par%owns_data = 0
    if (ag%w_par%size > 0) call set_cloc(ag%w_tmp, ag%w_par)

  end subroutine d_set_default_bcm_w

  subroutine set_cloc(vect,w_par)
    use iso_c_binding
    real(psb_dpk_), target :: vect(:)
    type(bcm_Vector) :: w_par
    
    w_par%data = c_loc(vect)
  end subroutine set_cloc


  subroutine d_bcmatch_bld_wnxt(ag,ilaggr,valaggr,nx)
    use psb_realloc_mod
    implicit none 
    class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_ipk_), intent(in)  :: ilaggr(:)
    real(psb_dpk_), intent(in)     :: valaggr(:)
    integer(psb_ipk_), intent(in)  :: nx

    integer(psb_ipk_) :: info,i,j

    call psb_realloc(nx,ag%w_nxt,info)
    associate(w_nxt => ag%w_nxt, w_tmp=>ag%w_tmp)
      w_nxt = dzero
      if (.false.) then 
        do j=1, size(ilaggr)
          i = ilaggr(j)
          w_nxt(i) = w_nxt(i) + valaggr(j)*w_tmp(j)
        end do
        !write(0,*) 'Old copy ',w_nxt(1:10)
      else
        !write(0,*) 'New copy ',nx
        do i=1, nx
          w_nxt(i) = w_tmp(i)
        end do
        !write(0,*) 'New copy ',w_nxt(1:10)
      end if
    end associate
    
  end subroutine d_bcmatch_bld_wnxt
  
  function d_bcmatch_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "BootCMatch aggregation"
  end function d_bcmatch_aggregator_fmt

  subroutine  d_bcmatch_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(mld_d_bcmatch_aggregator_type), intent(in) :: ag
    type(mld_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'BootCMatch Aggregator'
    write(iout,*) '   Number of BootCMatch sweeps: ',ag%n_sweeps
    write(iout,*) '   Matching algorithm         : ',ag%matching_alg
    write(iout,*) '    0: Preis 1: MC64  2: SPRAL  '
    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine d_bcmatch_aggregator_descr

  subroutine  d_bcmatch_aggregator_update_next(ag,agnext,info)
    use psb_realloc_mod
    implicit none 
    class(mld_d_bcmatch_aggregator_type), target, intent(inout) :: ag
    class(mld_d_base_aggregator_type), target, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    !
    !
    select type(agnext)
    class is (mld_d_bcmatch_aggregator_type)
      agnext%matching_alg = ag%matching_alg
      agnext%n_sweeps     = ag%n_sweeps 
      agnext%max_csize    = ag%max_csize 
      agnext%max_nlevels  = ag%max_nlevels
      ! Is this going to generate shallow copies/memory leaks/double frees?
      ! To be investigated further.
      call psb_safe_ab_cpy(ag%w_nxt,agnext%w_tmp,info)
      call agnext%set_c_default_w()
    class default
      ! What should we do here? 
    end select
    info = 0 
  end subroutine d_bcmatch_aggregator_update_next

  subroutine d_bcmatch_aggr_cseti(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(mld_d_bcmatch_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer, intent(in), optional                :: idx
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='d_bcmatch_aggr_cseti'
    info = psb_success_

    ! For now we ignore IDX
    
    select case(what)
    case('BCM_MATCH_ALG')
      ag%matching_alg=val
    case('BCM_SWEEPS')
      ag%n_sweeps=val
    case('BCM_MAX_CSIZE')
      ag%max_csize=val
    case('BCM_MAX_NLEVELS')
      ag%max_nlevels=val
    case('BCM_W_SIZE')
      call ag%bld_default_w(val)
    case default
      
    end select
    return
  end subroutine d_bcmatch_aggr_cseti

  subroutine d_bcmatch_aggr_set_default(ag)

    Implicit None

    ! Arguments
    class(mld_d_bcmatch_aggregator_type), intent(inout) :: ag
    character(len=20)  :: name='d_bcmatch_aggr_set_default'
    ag%matching_alg = 0
    ag%n_sweeps     = 1
    ag%max_nlevels  = 36
    ag%max_csize    = 10

    return

  end subroutine d_bcmatch_aggr_set_default

  subroutine  d_bcmatch_aggregator_free(ag,info)
    use iso_c_binding
    implicit none 
    class(mld_d_bcmatch_aggregator_type), intent(inout) :: ag
    integer(psb_ipk_), intent(out)       :: info

    info = 0
    if (allocated(ag%w_tmp)) deallocate(ag%w_tmp,stat=info)
    if (info /= 0) return
    if (allocated(ag%w_nxt)) deallocate(ag%w_nxt,stat=info)
    if (info /= 0) return
    ag%w_par%size = 0
    ag%w_par%data = c_null_ptr
    ag%w_par%owns_data = 0
  end subroutine d_bcmatch_aggregator_free

  subroutine  d_bcmatch_aggregator_clone(ag,agnext,info)
    implicit none 
    class(mld_d_bcmatch_aggregator_type), intent(inout) :: ag
    class(mld_d_base_aggregator_type), allocatable, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    info = 0 
    if (allocated(agnext)) then
      call agnext%free(info)
      if (info == 0) deallocate(agnext,stat=info)
    end if
    if (info /= 0) return
    allocate(agnext,source=ag,stat=info)
    select type(agnext)
    class is (mld_d_bcmatch_aggregator_type)
      call agnext%set_c_default_w()
    class default
      ! Should never ever get here
      info = -1
    end select
  end subroutine d_bcmatch_aggregator_clone  
  
end module mld_d_bcmatch_aggregator_mod

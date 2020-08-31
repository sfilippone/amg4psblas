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
!
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. 
!    
module amg_z_base_aggregator_mod

  use amg_base_prec_type, only : amg_dml_parms, amg_daggr_data 
  use psb_base_mod, only : psb_zspmat_type, psb_lzspmat_type, psb_z_vect_type, &
       & psb_z_base_vect_type, psb_zlinmap_type, psb_dpk_, &
       & psb_lz_csr_sparse_mat,  psb_lz_coo_sparse_mat, &
       & psb_z_csr_sparse_mat,  psb_z_coo_sparse_mat, &
       & psb_ipk_, psb_epk_, psb_lpk_, psb_desc_type, psb_i_base_vect_type, &
       & psb_erractionsave, psb_error_handler, psb_success_, psb_toupper
  !
  !  
  !
  !> \class amg_z_base_aggregator_type
  !!
  !!  It is the data type containing the basic interface definition for
  !!  building a multigrid hierarchy by aggregation. The base object has no attributes,
  !!  it is intended to be essentially an abstract type. 
  !!  
  !!
  !!  type amg_z_base_aggregator_type
  !!  end type
  !!  
  !!  
  !!  Methods:
  !!
  !!  bld_tprol   -   Build a tentative prolongator
  !!  
  !!  mat_bld     -   Build prolongator/restrictor and coarse matrix ac
  !!  
  !!  mat_asb     -   Convert  prolongator/restrictor/coarse matrix
  !!                  and fix their descriptor(s)
  !!                  
  !!  update_next -   Transfer information to the next level; default is
  !!                  to do nothing, i.e. aggregators at different
  !!                  levels are independent.
  !!  
  !!  default     -   Apply defaults
  !!  set_aggr_type - For aggregator that have internal options.
  !!  fmt         - Return a short string description
  !!  descr       - Print a more detailed description
  !!
  !!  cseti, csetr, csetc  - Set internal parameters, if any
  !  
  type amg_z_base_aggregator_type
    ! Do we want to purge explicit zeros when aggregating? 
    logical :: do_clean_zeros
  contains
    procedure, pass(ag) :: bld_tprol   => amg_z_base_aggregator_build_tprol
    procedure, pass(ag) :: mat_bld     => amg_z_base_aggregator_mat_bld
    procedure, pass(ag) :: mat_asb     => amg_z_base_aggregator_mat_asb
    procedure, pass(ag) :: bld_map     => amg_z_base_aggregator_bld_map
    procedure, pass(ag) :: update_next => amg_z_base_aggregator_update_next
    procedure, pass(ag) :: clone       => amg_z_base_aggregator_clone
    procedure, pass(ag) :: free        => amg_z_base_aggregator_free
    procedure, pass(ag) :: default     => amg_z_base_aggregator_default
    procedure, pass(ag) :: descr       => amg_z_base_aggregator_descr
    procedure, pass(ag) :: sizeof      => amg_z_base_aggregator_sizeof   
    procedure, pass(ag) :: set_aggr_type => amg_z_base_aggregator_set_aggr_type
    procedure, nopass   :: fmt         => amg_z_base_aggregator_fmt
    procedure, pass(ag) :: cseti       => amg_z_base_aggregator_cseti
    procedure, pass(ag) :: csetr       => amg_z_base_aggregator_csetr
    procedure, pass(ag) :: csetc       => amg_z_base_aggregator_csetc
    generic, public     :: set         => cseti, csetr, csetc
    procedure, nopass   :: xt_desc     => amg_z_base_aggregator_xt_desc
  end type amg_z_base_aggregator_type

  abstract interface  
    subroutine amg_z_soc_map_bld(iorder,theta,clean_zeros,a,desc_a,nlaggr,ilaggr,info)
      import :: psb_zspmat_type, psb_desc_type, psb_dpk_, psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_ipk_), intent(in)     :: iorder
      logical, intent(in)               :: clean_zeros
      type(psb_zspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_dpk_), intent(in)         :: theta
      integer(psb_lpk_), allocatable, intent(out)  :: ilaggr(:),nlaggr(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_z_soc_map_bld
  end interface

  interface amg_ptap
    subroutine amg_z_ptap(a_csr,desc_a,nlaggr,parms,ac,&
         & coo_prol,desc_cprol,coo_restr,info,desc_ax)
      import :: psb_z_csr_sparse_mat, psb_zspmat_type, psb_desc_type, &
           & psb_z_coo_sparse_mat, amg_dml_parms, psb_dpk_, psb_ipk_, psb_lpk_
      implicit none
      type(psb_z_csr_sparse_mat), intent(inout) :: a_csr
      type(psb_desc_type), intent(in)            :: desc_a
      integer(psb_lpk_), intent(inout)           :: nlaggr(:)
      type(amg_dml_parms), intent(inout)         :: parms 
      type(psb_z_coo_sparse_mat), intent(inout) :: coo_prol, coo_restr
      type(psb_desc_type), intent(inout)         :: desc_cprol
      type(psb_zspmat_type), intent(out)        :: ac
      integer(psb_ipk_), intent(out)             :: info
      type(psb_desc_type), intent(inout), optional :: desc_ax
    end subroutine amg_z_ptap
!!$    subroutine amg_z_lz_ptap(a_csr,desc_a,nlaggr,parms,ac,&
!!$         & coo_prol,desc_cprol,coo_restr,info,desc_ax)
!!$      import :: psb_z_csr_sparse_mat, psb_lzspmat_type, psb_desc_type, &
!!$           & psb_lz_coo_sparse_mat, amg_dml_parms, psb_dpk_, psb_ipk_, psb_lpk_
!!$      implicit none
!!$      type(psb_z_csr_sparse_mat), intent(inout) :: a_csr
!!$      type(psb_desc_type), intent(in)            :: desc_a
!!$      integer(psb_lpk_), intent(inout)           :: nlaggr(:)
!!$      type(amg_dml_parms), intent(inout)         :: parms 
!!$      type(psb_lz_coo_sparse_mat), intent(inout) :: coo_prol, coo_restr
!!$      type(psb_desc_type), intent(inout)         :: desc_cprol
!!$      type(psb_lzspmat_type), intent(out)        :: ac
!!$      integer(psb_ipk_), intent(out)             :: info
!!$      type(psb_desc_type), intent(inout), optional :: desc_ax
!!$    end subroutine amg_z_lz_ptap
!!$    subroutine amg_lz_ptap(a_csr,desc_a,nlaggr,parms,ac,&
!!$         & coo_prol,desc_cprol,coo_restr,info,desc_ax)
!!$      import :: psb_lz_csr_sparse_mat, psb_lzspmat_type, psb_desc_type, &
!!$           & psb_lz_coo_sparse_mat, amg_dml_parms, psb_dpk_, psb_ipk_, psb_lpk_
!!$      implicit none
!!$      type(psb_lz_csr_sparse_mat), intent(inout) :: a_csr
!!$      type(psb_desc_type), intent(in)            :: desc_a
!!$      integer(psb_lpk_), intent(inout)           :: nlaggr(:)
!!$      type(amg_dml_parms), intent(inout)         :: parms 
!!$      type(psb_lz_coo_sparse_mat), intent(inout) :: coo_prol, coo_restr
!!$      type(psb_desc_type), intent(inout)         :: desc_cprol
!!$      type(psb_lzspmat_type), intent(out)        :: ac
!!$      integer(psb_ipk_), intent(out)             :: info
!!$      type(psb_desc_type), intent(inout), optional :: desc_ax
!!$    end subroutine amg_lz_ptap
  end interface amg_ptap
  
contains

  subroutine amg_z_base_aggregator_cseti(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    ! Do nothing
    info = 0
  end subroutine amg_z_base_aggregator_cseti

  subroutine amg_z_base_aggregator_csetr(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    real(psb_dpk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    ! Do nothing
    info = 0
  end subroutine amg_z_base_aggregator_csetr

  subroutine amg_z_base_aggregator_csetc(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    character(len=*), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    ! Set clean zeros, or do nothing. 
    select case (psb_toupper(trim(what)))
    case('AGGR_CLEAN_ZEROS')
      select case (psb_toupper(trim(val)))
      case('TRUE','T')
        ag%do_clean_zeros = .true.
      case('FALSE','F')
        ag%do_clean_zeros = .false.
      end select
    end select
    info = 0
  end subroutine amg_z_base_aggregator_csetc

  
  subroutine  amg_z_base_aggregator_update_next(ag,agnext,info)
    implicit none 
    class(amg_z_base_aggregator_type), target, intent(inout) :: ag, agnext
    integer(psb_ipk_), intent(out)       :: info

    !
    ! Base version does nothing. 
    !
    info = 0 
  end subroutine amg_z_base_aggregator_update_next
  
  subroutine  amg_z_base_aggregator_clone(ag,agnext,info)
    implicit none 
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    class(amg_z_base_aggregator_type), allocatable, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    info = 0 
    if (allocated(agnext)) then
      call agnext%free(info)
      if (info == 0) deallocate(agnext,stat=info)
    end if
    if (info /= 0) return
    allocate(agnext,source=ag,stat=info)
    
  end subroutine amg_z_base_aggregator_clone

  subroutine  amg_z_base_aggregator_free(ag,info)
    implicit none 
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    integer(psb_ipk_), intent(out)       :: info
    
    info = psb_success_
    return
  end subroutine amg_z_base_aggregator_free
  
  subroutine  amg_z_base_aggregator_default(ag)
    implicit none 
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    ! Only one default setting
    ag%do_clean_zeros = .true.
    
    return
  end subroutine amg_z_base_aggregator_default

  function amg_z_base_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "Default aggregator "
  end function amg_z_base_aggregator_fmt

  function amg_z_base_aggregator_sizeof(ag) result(val)
    implicit none
    class(amg_z_base_aggregator_type), intent(in)  :: ag
    integer(psb_epk_)  :: val

    val = 1
  end function amg_z_base_aggregator_sizeof

  function amg_z_base_aggregator_xt_desc() result(val)
    implicit none 
    logical  :: val

    val = .false.
  end function amg_z_base_aggregator_xt_desc

  subroutine  amg_z_base_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(amg_z_base_aggregator_type), intent(in) :: ag
    type(amg_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine amg_z_base_aggregator_descr
  
  subroutine  amg_z_base_aggregator_set_aggr_type(ag,parms,info)
    implicit none 
    class(amg_z_base_aggregator_type), intent(inout) :: ag
    type(amg_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(out) :: info

    ! Do nothing
    
    return
  end subroutine amg_z_base_aggregator_set_aggr_type

  !
  !> Function  bld_tprol:
  !! \memberof  amg_z_base_aggregator_type
  !! \brief  Build a tentative prolongator.           
  !!         The routine will map the local matrix entries to aggregates.
  !!         The mapping is store in ILAGGR; for each local row index I,      
  !!         ILAGGR(I) contains the index of the aggregate to which index I
  !!         will contribute, in global numbering. 
  !!         Many aggregations produce a binary tentative prolongator, but some
  !!         do not, hence we also need the OP_PROL output.
  !!         AG_DATA is passed here just in case some of the
  !!         aggregators need it internally, most of them will ignore.
  !!         
  !!  \param ag      The input aggregator object
  !!  \param parms   The auxiliary parameters object
  !!  \param ag_data Auxiliary global aggregation info
  !!  \param a       The local matrix part
  !!  \param desc_a  The descriptor
  !!  \param ilaggr  Output aggregation map
  !!  \param nlaggr  Sizes of ilaggr on all processes
  !!  \param op_prol The tentative prolongator operator
  !!  \param info    Return code
  !!           
  !
  subroutine  amg_z_base_aggregator_build_tprol(ag,parms,ag_data,&
       & a,desc_a,ilaggr,nlaggr,t_prol,info)
    use psb_base_mod
    implicit none
    class(amg_z_base_aggregator_type), target, intent(inout) :: ag
    type(amg_dml_parms), intent(inout)  :: parms
    type(amg_daggr_data), intent(in)    :: ag_data
    type(psb_zspmat_type), intent(inout) :: a
    type(psb_desc_type), intent(inout)     :: desc_a
    integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
    type(psb_lzspmat_type), intent(out)  :: t_prol
    integer(psb_ipk_), intent(out)      :: info
    
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='z_base_aggregator_build_tprol'
    
    call psb_erractionsave(err_act)
    
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return
    
9999 call psb_error_handler(err_act)
    
    return

  end subroutine amg_z_base_aggregator_build_tprol

  !
  !> Function   mat_bld    
  !! \memberof  amg_z_base_aggregator_type
  !! \brief     Build prolongator/restrictor/coarse matrix.
  !!
  !!
  !!  \param ag        The input aggregator object
  !!  \param parms     The auxiliary parameters object
  !!  \param a         The local matrix part
  !!  \param desc_a    The descriptor
  !!  \param ilaggr    Aggregation map
  !!  \param nlaggr    Sizes of ilaggr on all processes
  !!  \param ac        On output the coarse matrix
  !!  \param op_prol   On input, the  tentative prolongator operator, on output
  !!                   the final prolongator
  !!  \param op_restr  On output, the restrictor operator;
  !!                   in many cases it is the transpose of the prolongator. 
  !!  \param info    Return code
  !!  
  subroutine  amg_z_base_aggregator_mat_bld(ag,parms,a,desc_a,ilaggr,nlaggr,&
       & ac,desc_ac,op_prol,op_restr,t_prol,info)
    use psb_base_mod
    implicit none
    class(amg_z_base_aggregator_type), target, intent(inout) :: ag
    type(amg_dml_parms), intent(inout)   :: parms 
    type(psb_zspmat_type), intent(in)    :: a
    type(psb_desc_type), intent(inout)   :: desc_a
    integer(psb_lpk_), intent(inout)     :: ilaggr(:), nlaggr(:)
    type(psb_lzspmat_type), intent(inout) :: t_prol
    type(psb_zspmat_type), intent(out)   :: op_prol, ac,op_restr
    type(psb_desc_type), intent(inout)   :: desc_ac
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='z_base_aggregator_mat_bld'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine amg_z_base_aggregator_mat_bld

  !
  !> Function   mat_asb    
  !! \memberof  amg_z_base_aggregator_type
  !! \brief     Build prolongator/restrictor/coarse matrix.
  !!
  !!
  !!  \param ag        The input aggregator object
  !!  \param parms     The auxiliary parameters object
  !!  \param a         The local matrix part
  !!  \param desc_a    The descriptor
  !!  \param ilaggr    Aggregation map
  !!  \param nlaggr    Sizes of ilaggr on all processes
  !!  \param ac        On output the coarse matrix
  !!  \param op_prol   On input, the  tentative prolongator operator, on output
  !!                   the final prolongator
  !!  \param op_restr  On output, the restrictor operator;
  !!                   in many cases it is the transpose of the prolongator. 
  !!  \param info    Return code
  !!  
  subroutine  amg_z_base_aggregator_mat_asb(ag,parms,a,desc_a,&
       & ac,desc_ac, op_prol,op_restr,info)
    use psb_base_mod
    implicit none
    class(amg_z_base_aggregator_type), target, intent(inout) :: ag
    type(amg_dml_parms), intent(inout)   :: parms 
    type(psb_zspmat_type), intent(in)    :: a
    type(psb_desc_type), intent(inout)     :: desc_a
    type(psb_zspmat_type), intent(inout) :: op_prol,ac,op_restr
    type(psb_desc_type), intent(inout)      :: desc_ac
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='z_base_aggregator_mat_asb'

    call psb_erractionsave(err_act)

    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine amg_z_base_aggregator_mat_asb

  !
  !> Function   bld_map
  !! \memberof  amg_z_base_aggregator_type
  !! \brief     Build linear map between hierarchy levels 
  !!
  !!
  !!  \param ag        The input aggregator object
  !!  \param desc_a    The fine space descriptor
  !!  \param desc_ac   The coarse space descriptor
  !!  \param ilaggr    Aggregation map vector
  !!  \param nlaggr    Sizes of ilaggr on all processes
  !!  \param op_prol   The  prolongator operator
  !!  \param op_restr  The restrictor operator
  !!  \param map       The output map
  !!  \param info    Return code
  !!  
  subroutine  amg_z_base_aggregator_bld_map(ag,desc_a,desc_ac,ilaggr,nlaggr,&
       & op_restr,op_prol,map,info)
    use psb_base_mod
    implicit none
    class(amg_z_base_aggregator_type), target, intent(inout) :: ag
    type(psb_desc_type), intent(in), target :: desc_a, desc_ac
    integer(psb_lpk_), intent(inout)     :: ilaggr(:), nlaggr(:)
    type(psb_zspmat_type), intent(inout)   :: op_restr, op_prol
    type(psb_zlinmap_type), intent(out)    :: map
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='z_base_aggregator_bld_map'

    call psb_erractionsave(err_act)
    !
    ! Copy the prolongation/restriction matrices into the descriptor map.
    !  op_restr => PR^T   i.e. restriction  operator
    !  op_prol => PR     i.e. prolongation operator
    !
    !  WARNING: need to check whether the  copy into IOP_RESTR/IOP_PROL
    !  is safe or not.
    !
    !  This default implementation reuses desc_a/desc_ac through
    !  pointers in the map structure.
    !      
    map = psb_linmap(psb_map_aggr_,desc_a,&
         & desc_ac,op_restr,op_prol,ilaggr,nlaggr)

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine amg_z_base_aggregator_bld_map


end module amg_z_base_aggregator_mod

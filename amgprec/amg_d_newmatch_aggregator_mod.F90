!
!
!  The aggregator object hosts the aggregation method for building
!  the multilevel hierarchy. This variant is based on the hybrid method
!  presented in 
!
!
!   sm           -  class(amg_T_base_smoother_type), allocatable
!                   The current level preconditioner (aka smoother).
!   parms        -  type(amg_RTml_parms)
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

module amg_d_newmatch_aggregator_mod
  use amg_d_base_aggregator_mod
  use iso_c_binding
  
  type, bind(c)::  nwm_Vector
    type(c_ptr) :: data
    integer(c_int) :: size
    integer(c_int) :: owns_data
  end type nwm_Vector
  
  type, bind(c)::  nwm_CSRMatrix
    type(c_ptr) :: i
    type(c_ptr) :: j
    integer(c_int) :: num_rows
    integer(c_int) :: num_cols
    integer(c_int) :: num_nonzeros
    integer(c_int) :: owns_data
    type(c_ptr) :: data
  end type nwm_CSRMatrix
  
  type, extends(amg_d_base_aggregator_type) :: amg_d_newmatch_aggregator_type
    integer(psb_ipk_) :: matching_alg
    integer(psb_ipk_) :: n_sweeps
    !
    !  Note: the BootCMatch kernel we invoke  overwrites
    !  the W argument with its update. Hence, copy it in w_nxt
    !  before passing it to the matching
    !
    integer(psb_ipk_) :: orig_aggr_size
    integer(psb_ipk_) :: jacobi_sweeps
    real(psb_dpk_), allocatable :: w(:), w_nxt(:)
    type(psb_dspmat_type), allocatable  :: prol, restr
    type(psb_dspmat_type), allocatable  :: ac, base_a, rwa
    type(psb_desc_type), allocatable    :: desc_ac, desc_ax, base_desc, rwdesc
    type(nwm_Vector)  :: w_c_nxt
    integer(psb_ipk_) :: max_csize
    integer(psb_ipk_) :: max_nlevels
    logical           :: reproducible_matching = .false.
    logical           :: need_symmetrize       = .false.
    logical           :: unsmoothed_hierarchy  = .true.
  contains
    procedure, pass(ag) :: bld_tprol    => amg_d_newmatch_aggregator_build_tprol
    procedure, pass(ag) :: cseti        => d_newmatch_aggr_cseti
    procedure, pass(ag) :: default      => d_newmatch_aggr_set_default
    procedure, pass(ag) :: mat_asb      => amg_d_newmatch_aggregator_mat_asb
    procedure, pass(ag) :: mat_bld      => amg_d_newmatch_aggregator_mat_bld
    procedure, pass(ag) :: inner_mat_asb => amg_d_newmatch_aggregator_inner_mat_asb
    procedure, pass(ag) :: update_next  => d_newmatch_aggregator_update_next
    procedure, pass(ag) :: bld_wnxt     => d_newmatch_bld_wnxt
    procedure, pass(ag) :: bld_default_w    => d_bld_default_w
    procedure, pass(ag) :: set_c_default_w  => d_set_default_nwm_w
    procedure, pass(ag) :: descr        => d_newmatch_aggregator_descr
    procedure, pass(ag) :: clone        => d_newmatch_aggregator_clone
    procedure, pass(ag) :: free         => d_newmatch_aggregator_free
    procedure, nopass   :: fmt          => d_newmatch_aggregator_fmt
  end type amg_d_newmatch_aggregator_type


  interface
    subroutine  amg_d_newmatch_aggregator_build_tprol(ag,parms,ag_data,&
         & a,desc_a,ilaggr,nlaggr,t_prol,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, &
           & psb_dspmat_type, psb_ldspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_lpk_, psb_epk_, amg_dml_parms, amg_daggr_data
      implicit none
      class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_dml_parms), intent(inout)   :: parms 
      type(amg_daggr_data), intent(in)     :: ag_data
      type(psb_dspmat_type), intent(inout) :: a
      type(psb_desc_type), intent(inout)   :: desc_a
      integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_ldspmat_type), intent(out)  :: t_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine amg_d_newmatch_aggregator_build_tprol
  end interface

  interface
    subroutine  amg_d_newmatch_aggregator_mat_bld(ag,parms,a,desc_a,ilaggr,nlaggr,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, &
           & psb_dspmat_type, psb_ldspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_lpk_, psb_epk_, amg_dml_parms
      implicit none
      class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_dml_parms), intent(inout)    :: parms 
      type(psb_dspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(inout)    :: desc_a
      integer(psb_lpk_), intent(inout)      :: ilaggr(:), nlaggr(:)
      type(psb_dspmat_type), intent(out)    :: op_prol,ac,op_restr
      type(psb_ldspmat_type), intent(inout) :: t_prol
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)        :: info
    end subroutine amg_d_newmatch_aggregator_mat_bld
  end interface

  interface
    subroutine  amg_d_newmatch_aggregator_mat_asb(ag,parms,a,desc_a,&
         & ac,desc_ac, op_prol,op_restr,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, &
           & psb_dspmat_type, psb_ldspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_lpk_, psb_epk_, amg_dml_parms
      implicit none
      class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_dml_parms), intent(inout)    :: parms 
      type(psb_dspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(inout)    :: desc_a
      type(psb_dspmat_type), intent(inout) :: op_prol,ac,op_restr
      type(psb_desc_type), intent(inout)     :: desc_ac
      integer(psb_ipk_), intent(out)        :: info
    end subroutine amg_d_newmatch_aggregator_mat_asb
  end interface  
  

  interface
    subroutine amg_d_newmatch_map_to_tprol(desc_a,ilaggr,nlaggr,valaggr, op_prol,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, &
           & psb_dspmat_type, psb_ldspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_lpk_, psb_epk_, amg_dml_parms
      implicit none
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_lpk_), allocatable, intent(inout)  :: ilaggr(:),nlaggr(:)
      real(psb_dpk_), allocatable, intent(inout)  :: valaggr(:)
      type(psb_ldspmat_type), intent(out)  :: op_prol
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_d_newmatch_map_to_tprol
  end interface
  
  interface
    subroutine amg_daggrmat_unsmth_spmm_asb(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,op_prol,op_restr,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, &
           & psb_dspmat_type, psb_ldspmat_type, psb_dpk_,  &
           & psb_ipk_, psb_lpk_, psb_epk_, amg_dml_parms
      implicit none
      type(psb_dspmat_type), intent(in)        :: a
      type(psb_desc_type), intent(in)            :: desc_a
      integer(psb_lpk_), intent(inout)           :: ilaggr(:), nlaggr(:)
      type(amg_dml_parms), intent(inout)      :: parms 
      type(psb_ldspmat_type), intent(inout)     :: op_prol
      type(psb_ldspmat_type), intent(out)       :: ac,op_restr
      integer(psb_ipk_), intent(out)             :: info
    end subroutine amg_daggrmat_unsmth_spmm_asb
  end interface

  interface
    subroutine  amg_d_newmatch_aggregator_inner_mat_asb(ag,parms,a,desc_a,&
         & ac,desc_ac, op_prol,op_restr,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, psb_dspmat_type,&
           & psb_ldspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_dml_parms, amg_daggr_data
      implicit none
      class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_dml_parms), intent(inout)    :: parms
      type(psb_dspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(in)       :: desc_a
      type(psb_dspmat_type), intent(inout) :: op_prol,op_restr
      type(psb_dspmat_type), intent(inout)  :: ac
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)        :: info
    end subroutine amg_d_newmatch_aggregator_inner_mat_asb
  end interface

  
!!$  interface
!!$    subroutine amg_d_newmatch_unsmth_spmm_bld(a,desc_a,ilaggr,nlaggr,parms,&
!!$         & ac,desc_ac,op_prol,op_restr,t_prol,info)
!!$      import :: amg_d_newmatch_aggregator_type, psb_desc_type, &
!!$           & psb_dspmat_type, psb_ldspmat_type, psb_dpk_,  &
!!$           & psb_ipk_, psb_lpk_, psb_epk_, amg_dml_parms
!!$      
!!$      implicit none
!!$      
!!$      ! Arguments
!!$      type(psb_dspmat_type), intent(in)        :: a
!!$      type(psb_desc_type), intent(in)            :: desc_a
!!$      integer(psb_lpk_), intent(inout)           :: ilaggr(:), nlaggr(:)
!!$      type(amg_dml_parms), intent(inout)      :: parms 
!!$      type(psb_ldspmat_type), intent(inout)     :: t_prol
!!$      type(psb_dspmat_type), intent(inout)       :: op_prol,ac,op_restr
!!$      type(psb_desc_type), intent(inout)    :: desc_ac
!!$      integer(psb_ipk_), intent(out)             :: info
!!$    end subroutine amg_d_newmatch_unsmth_spmm_bld
!!$  end interface

  interface
    subroutine amg_d_newmatch_spmm_bld_ov(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, psb_dspmat_type,&
           & psb_ldspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_dml_parms, amg_daggr_data
      implicit none
      type(psb_dspmat_type), intent(inout)    :: a
      type(psb_desc_type), intent(inout)       :: desc_a
      integer(psb_lpk_), intent(inout)        :: ilaggr(:), nlaggr(:)
      type(amg_dml_parms), intent(inout)      :: parms
      type(psb_ldspmat_type), intent(inout)   :: t_prol
      type(psb_dspmat_type), intent(inout)    :: op_prol,ac, op_restr
      type(psb_desc_type), intent(out)          :: desc_ac
      integer(psb_ipk_), intent(out)          :: info
    end subroutine amg_d_newmatch_spmm_bld_ov
  end interface

  interface
    subroutine amg_d_newmatch_spmm_bld_inner(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_d_newmatch_aggregator_type, psb_desc_type, psb_dspmat_type,&
           & psb_ldspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_dml_parms, amg_daggr_data,&
           & psb_d_csr_sparse_mat, psb_ld_csr_sparse_mat
      implicit none
      type(psb_d_csr_sparse_mat), intent(inout) :: a
      type(psb_desc_type), intent(inout)          :: desc_a
      integer(psb_lpk_), intent(inout)          :: ilaggr(:), nlaggr(:)
      type(amg_dml_parms), intent(inout)        :: parms
      type(psb_ldspmat_type), intent(inout)     :: t_prol
      type(psb_dspmat_type), intent(inout)      :: op_prol,ac, op_restr
      type(psb_desc_type), intent(out)          :: desc_ac
      integer(psb_ipk_), intent(out)            :: info
    end subroutine amg_d_newmatch_spmm_bld_inner
  end interface
  
  private :: is_legal_malg, is_legal_csize, is_legal_nsweeps, is_legal_nlevels
  

contains

  subroutine d_bld_default_w(ag,nr)
    use psb_realloc_mod
    implicit none 
    class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_ipk_), intent(in) :: nr
    integer(psb_ipk_) :: info
    call psb_realloc(nr,ag%w,info)
    if (info /= psb_success_) return
    ag%w = done
    call ag%set_c_default_w()
  end subroutine d_bld_default_w

  subroutine d_set_default_nwm_w(ag)
    use psb_realloc_mod
    use iso_c_binding
    implicit none 
    class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_ipk_) :: info
    
    call psb_safe_ab_cpy(ag%w,ag%w_nxt,info)
    ag%w_c_nxt%size      = psb_size(ag%w_nxt)
    ag%w_c_nxt%owns_data = 0
    if (ag%w_c_nxt%size > 0) call set_cloc(ag%w_nxt, ag%w_c_nxt)

  end subroutine d_set_default_nwm_w

  subroutine set_cloc(vect,w_c_nxt)
    use iso_c_binding
    real(psb_dpk_), target :: vect(:)
    type(nwm_Vector) :: w_c_nxt
    
    w_c_nxt%data = c_loc(vect)
  end subroutine set_cloc


  subroutine d_newmatch_bld_wnxt(ag,ilaggr,valaggr,nx)
    use psb_realloc_mod
    implicit none 
    class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_lpk_), intent(in)  :: ilaggr(:)
    real(psb_dpk_), intent(in)     :: valaggr(:)
    integer(psb_ipk_), intent(in)  :: nx

    integer(psb_ipk_) :: info,i,j

    ! The vector was already fixed in the call to Newmatch.
    call psb_realloc(nx,ag%w_nxt,info)
    
  end subroutine d_newmatch_bld_wnxt
  
  function d_newmatch_aggregator_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "new matching aggregation"
  end function d_newmatch_aggregator_fmt

  subroutine  d_newmatch_aggregator_descr(ag,parms,iout,info)
    implicit none 
    class(amg_d_newmatch_aggregator_type), intent(in) :: ag
    type(amg_dml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'NewMatch Aggregator'
    write(iout,*) '   Number of Matching   sweeps: ',ag%n_sweeps
    write(iout,*) '   Matching algorithm         : ',ag%matching_alg
    !write(iout,*) '    0: Preis 1: MC64  2: SPRAL  '
    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)
    
    return
  end subroutine d_newmatch_aggregator_descr

  function is_legal_malg(alg) result(val)
    logical :: val
    integer(psb_ipk_) :: alg

    val = ((0<=alg).and.(alg<=2))
  end function is_legal_malg

  function is_legal_csize(csize) result(val)
    logical :: val
    integer(psb_ipk_) :: csize

    val = ((-1==csize).or.(csize >0))
  end function is_legal_csize

  function is_legal_nsweeps(nsw) result(val)
    logical :: val
    integer(psb_ipk_) :: nsw

    val = (1<=nsw)
  end function is_legal_nsweeps

  function is_legal_nlevels(nlv) result(val)
    logical :: val
    integer(psb_ipk_) :: nlv

    val = (1<=nlv)
  end function is_legal_nlevels
  
  subroutine  d_newmatch_aggregator_update_next(ag,agnext,info)
    use psb_realloc_mod
    implicit none 
    class(amg_d_newmatch_aggregator_type), target, intent(inout) :: ag
    class(amg_d_base_aggregator_type), target, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    !
    !
    select type(agnext)
    class is (amg_d_newmatch_aggregator_type)
      if (.not.is_legal_malg(agnext%matching_alg)) &
           & agnext%matching_alg = ag%matching_alg
      if (.not.is_legal_nsweeps(agnext%n_sweeps))&
           & agnext%n_sweeps     = ag%n_sweeps 
      if (.not.is_legal_csize(agnext%max_csize))&
           & agnext%max_csize    = ag%max_csize 
      if (.not.is_legal_nlevels(agnext%max_nlevels))&
           & agnext%max_nlevels  = ag%max_nlevels
      ! Is this going to generate shallow copies/memory leaks/double frees?
      ! To be investigated further.
      call psb_safe_ab_cpy(ag%w_nxt,agnext%w,info)
      call agnext%set_c_default_w()
    class default
      ! What should we do here? 
    end select
    info = 0 
  end subroutine d_newmatch_aggregator_update_next

  subroutine d_newmatch_aggr_cseti(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_newmatch_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='d_newmatch_aggr_cseti'
    info = psb_success_

    ! For now we ignore IDX
    
    select case(psb_toupper(trim(what)))
    case('NWM_MATCH_ALG')
      ag%matching_alg=val
    case('NWM_SWEEPS')
      ag%n_sweeps=val
    case('NWM_MAX_CSIZE')
      ag%max_csize=val
    case('NWM_MAX_NLEVELS')
      ag%max_nlevels=val
    case('NWM_W_SIZE')
      call ag%bld_default_w(val)
    case('AGGR_SIZE')
      ag%orig_aggr_size = val
      ag%n_sweeps=max(1,ceiling(log(val*1.0)/log(2.0)))
    case default
      ! Do nothing
    end select
    return
  end subroutine d_newmatch_aggr_cseti

  subroutine d_newmatch_aggr_set_default(ag)

    Implicit None

    ! Arguments
    class(amg_d_newmatch_aggregator_type), intent(inout) :: ag
    character(len=20)  :: name='d_newmatch_aggr_set_default'
    ag%matching_alg = 0
    ag%n_sweeps     = 1
    ag%max_nlevels  = 36
    ag%max_csize    = -1
    !
    ! Apparently newMatch works better
    ! by keeping all entries
    ! 
    ag%do_clean_zeros = .false.

    return

  end subroutine d_newmatch_aggr_set_default

  subroutine  d_newmatch_aggregator_free(ag,info)
    use iso_c_binding
    implicit none 
    class(amg_d_newmatch_aggregator_type), intent(inout) :: ag
    integer(psb_ipk_), intent(out)       :: info

    info = 0
    if (allocated(ag%w)) deallocate(ag%w,stat=info)
    if (info /= 0) return
    if (allocated(ag%w_nxt)) deallocate(ag%w_nxt,stat=info)
    if (info /= 0) return
    ag%w_c_nxt%size = 0
    ag%w_c_nxt%data = c_null_ptr
    ag%w_c_nxt%owns_data = 0
  end subroutine d_newmatch_aggregator_free

  subroutine  d_newmatch_aggregator_clone(ag,agnext,info)
    implicit none 
    class(amg_d_newmatch_aggregator_type), intent(inout) :: ag
    class(amg_d_base_aggregator_type), allocatable, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    info = 0 
    if (allocated(agnext)) then
      call agnext%free(info)
      if (info == 0) deallocate(agnext,stat=info)
    end if
    if (info /= 0) return
    allocate(agnext,source=ag,stat=info)
    select type(agnext)
    class is (amg_d_newmatch_aggregator_type)
      call agnext%set_c_default_w()
    class default
      ! Should never ever get here
      info = -1
    end select
  end subroutine d_newmatch_aggregator_clone  
  
end module amg_d_newmatch_aggregator_mod

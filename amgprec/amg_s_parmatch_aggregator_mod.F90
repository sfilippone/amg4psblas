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
! moved here from amg4psblas-extension
!
!
!                             AMG4PSBLAS  Extensions
!
!    (C) Copyright 2019
!
!                        Salvatore Filippone  Cranfield University
!        Pasqua D'Ambra         IAC-CNR, Naples, IT
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

module amg_s_parmatch_aggregator_mod
  use amg_s_base_aggregator_mod
  use smatchboxp_mod
#if defined(SERIAL_MPI)
  type, extends(amg_s_base_aggregator_type) :: amg_s_parmatch_aggregator_type
  end type amg_s_parmatch_aggregator_type
#else 
  type, extends(amg_s_base_aggregator_type) :: amg_s_parmatch_aggregator_type
    integer(psb_ipk_) :: matching_alg
    integer(psb_ipk_) :: n_sweeps   ! When n_sweeps >1 we need an auxiliary descriptor
    integer(psb_ipk_) :: orig_aggr_size
    integer(psb_ipk_) :: jacobi_sweeps
    real(psb_spk_), allocatable :: w(:), w_nxt(:)
    type(psb_sspmat_type), allocatable  :: prol, restr
    type(psb_sspmat_type), allocatable  :: ac, base_a, rwa
    type(psb_desc_type), allocatable    :: desc_ac, desc_ax, base_desc, rwdesc
    integer(psb_ipk_) :: max_csize
    integer(psb_ipk_) :: max_nlevels
    logical           :: reproducible_matching = .false.
    logical           :: need_symmetrize       = .false.
    logical           :: unsmoothed_hierarchy  = .true.
  contains
    procedure, pass(ag) :: bld_tprol    => amg_s_parmatch_aggregator_build_tprol
    procedure, pass(ag) :: mat_bld      => amg_s_parmatch_aggregator_mat_bld
    procedure, pass(ag) :: mat_asb      => amg_s_parmatch_aggregator_mat_asb
    procedure, pass(ag) :: inner_mat_asb      => amg_s_parmatch_aggregator_inner_mat_asb
    procedure, pass(ag) :: bld_map      => amg_s_parmatch_aggregator_bld_map
    procedure, pass(ag) :: csetc        => s_parmatch_aggr_csetc
    procedure, pass(ag) :: cseti        => s_parmatch_aggr_cseti
    procedure, pass(ag) :: default      => s_parmatch_aggr_set_default
    procedure, pass(ag) :: sizeof       => s_parmatch_aggregator_sizeof
    procedure, pass(ag) :: update_next  => s_parmatch_aggregator_update_next
    procedure, pass(ag) :: bld_wnxt     => s_parmatch_bld_wnxt
    procedure, pass(ag) :: bld_default_w    => s_bld_default_w
    procedure, pass(ag) :: set_c_default_w  => s_set_prm_c_default_w
    procedure, pass(ag) :: descr        => s_parmatch_aggregator_descr
    procedure, pass(ag) :: clone        => s_parmatch_aggregator_clone
    procedure, pass(ag) :: free         => s_parmatch_aggregator_free
    procedure, nopass   :: fmt          => s_parmatch_aggregator_fmt
    procedure, nopass   :: xt_desc      => amg_s_parmatch_aggregator_xt_desc
  end type amg_s_parmatch_aggregator_type

  interface
    subroutine  amg_s_parmatch_aggregator_build_tprol(ag,parms,ag_data,&
         & a,desc_a,ilaggr,nlaggr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_sml_parms), intent(inout)   :: parms
      type(amg_saggr_data), intent(in)     :: ag_data
      type(psb_sspmat_type), intent(inout) :: a
      type(psb_desc_type), intent(inout)     :: desc_a
      integer(psb_lpk_), allocatable, intent(out) :: ilaggr(:),nlaggr(:)
      type(psb_lsspmat_type), intent(out)  :: t_prol
      integer(psb_ipk_), intent(out)      :: info
    end subroutine amg_s_parmatch_aggregator_build_tprol
  end interface

  interface
    subroutine  amg_s_parmatch_aggregator_mat_bld(ag,parms,a,desc_a,ilaggr,nlaggr,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_sml_parms), intent(inout)   :: parms
      type(psb_sspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(inout)   :: desc_a
      integer(psb_lpk_), intent(inout)     :: ilaggr(:), nlaggr(:)
      type(psb_lsspmat_type), intent(inout)   :: t_prol
      type(psb_sspmat_type), intent(out)   ::  op_prol,ac,op_restr
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)       :: info
    end subroutine amg_s_parmatch_aggregator_mat_bld
  end interface

  interface
    subroutine  amg_s_parmatch_aggregator_mat_asb(ag,parms,a,desc_a,&
         & ac,desc_ac, op_prol,op_restr,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_sml_parms), intent(inout)    :: parms
      type(psb_sspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(inout)    :: desc_a
      type(psb_sspmat_type), intent(inout) :: op_prol,ac,op_restr
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)        :: info
    end subroutine amg_s_parmatch_aggregator_mat_asb
  end interface


  interface
    subroutine  amg_s_parmatch_aggregator_inner_mat_asb(ag,parms,a,desc_a,&
         & ac,desc_ac, op_prol,op_restr,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
      type(amg_sml_parms), intent(inout)    :: parms
      type(psb_sspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(in)       :: desc_a
      type(psb_sspmat_type), intent(inout) :: op_prol,op_restr
      type(psb_sspmat_type), intent(inout)  :: ac
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)        :: info
    end subroutine amg_s_parmatch_aggregator_inner_mat_asb
  end interface


  interface
    subroutine amg_s_parmatch_spmm_bld(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      type(psb_sspmat_type), intent(in)       :: a
      type(psb_desc_type), intent(inout)       :: desc_a
      integer(psb_lpk_), intent(inout)        :: ilaggr(:), nlaggr(:)
      type(amg_sml_parms), intent(inout)      :: parms
      type(psb_lsspmat_type), intent(inout)   :: t_prol
      type(psb_sspmat_type), intent(inout)    :: ac, op_prol, op_restr
      type(psb_desc_type), intent(out)        :: desc_ac
      integer(psb_ipk_), intent(out)          :: info
    end subroutine amg_s_parmatch_spmm_bld
  end interface

  interface
    subroutine amg_s_parmatch_unsmth_bld(ag,a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
      type(psb_sspmat_type), intent(in)       :: a
      type(psb_desc_type), intent(inout)      :: desc_a
      integer(psb_lpk_), intent(inout)        :: ilaggr(:), nlaggr(:)
      type(amg_sml_parms), intent(inout)      :: parms
      type(psb_lsspmat_type), intent(inout)   :: t_prol
      type(psb_sspmat_type), intent(inout)      :: op_prol,ac, op_restr
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)          :: info
    end subroutine amg_s_parmatch_unsmth_bld
  end interface

  interface
    subroutine amg_s_parmatch_smth_bld(ag,a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
      type(psb_sspmat_type), intent(in)       :: a
      type(psb_desc_type), intent(inout)      :: desc_a
      integer(psb_lpk_), intent(inout)        :: ilaggr(:), nlaggr(:)
      type(amg_sml_parms), intent(inout)      :: parms
      type(psb_lsspmat_type), intent(inout)   :: t_prol
      type(psb_sspmat_type), intent(inout)    :: op_prol,ac, op_restr
      type(psb_desc_type), intent(inout)    :: desc_ac
      integer(psb_ipk_), intent(out)          :: info
    end subroutine amg_s_parmatch_smth_bld
  end interface

  interface
    subroutine amg_s_parmatch_spmm_bld_ov(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data
      implicit none
      type(psb_sspmat_type), intent(inout)    :: a
      type(psb_desc_type), intent(inout)       :: desc_a
      integer(psb_lpk_), intent(inout)        :: ilaggr(:), nlaggr(:)
      type(amg_sml_parms), intent(inout)      :: parms
      type(psb_lsspmat_type), intent(inout)   :: t_prol
      type(psb_sspmat_type), intent(inout)    :: op_prol,ac, op_restr
      type(psb_desc_type), intent(out)          :: desc_ac
      integer(psb_ipk_), intent(out)          :: info
    end subroutine amg_s_parmatch_spmm_bld_ov
  end interface

  interface
    subroutine amg_s_parmatch_spmm_bld_inner(a,desc_a,ilaggr,nlaggr,parms,&
         & ac,desc_ac,op_prol,op_restr,t_prol,info)
      import :: amg_s_parmatch_aggregator_type, psb_desc_type, psb_sspmat_type,&
           & psb_lsspmat_type, psb_dpk_, psb_ipk_, psb_lpk_, amg_sml_parms, amg_saggr_data,&
           & psb_s_csr_sparse_mat, psb_ls_csr_sparse_mat
      implicit none
      type(psb_s_csr_sparse_mat), intent(inout) :: a
      type(psb_desc_type), intent(inout)          :: desc_a
      integer(psb_lpk_), intent(inout)          :: ilaggr(:), nlaggr(:)
      type(amg_sml_parms), intent(inout)        :: parms
      type(psb_lsspmat_type), intent(inout)     :: t_prol
      type(psb_sspmat_type), intent(inout)      :: op_prol,ac, op_restr
      type(psb_desc_type), intent(out)          :: desc_ac
      integer(psb_ipk_), intent(out)            :: info
    end subroutine amg_s_parmatch_spmm_bld_inner
  end interface

  private :: is_legal_malg, is_legal_csize, is_legal_nsweeps, is_legal_nlevels

contains

  subroutine s_bld_default_w(ag,nr)
    use psb_realloc_mod
    implicit none
    class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_ipk_), intent(in) :: nr
    integer(psb_ipk_) :: info
    call psb_realloc(nr,ag%w,info)
    if (info /= psb_success_) return
    ag%w = done
    !call ag%set_c_default_w()
  end subroutine s_bld_default_w

  subroutine s_set_prm_c_default_w(ag)
    use psb_realloc_mod
    use iso_c_binding
    implicit none
    class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_ipk_) :: info

    !write(0,*) 'prm_c_deafult_w '
    call psb_safe_ab_cpy(ag%w,ag%w_nxt,info)

  end subroutine s_set_prm_c_default_w

  subroutine s_parmatch_bld_wnxt(ag,ilaggr,valaggr,nx)
    use psb_realloc_mod
    implicit none
    class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
    integer(psb_lpk_), intent(in)  :: ilaggr(:)
    real(psb_spk_), intent(in)     :: valaggr(:)
    integer(psb_ipk_), intent(in)  :: nx

    integer(psb_ipk_) :: info,i,j

    ! The vector was already fixed in the call to BCMatch.
    !write(0,*) 'Executing bld_wnxt ',nx
    call psb_realloc(nx,ag%w_nxt,info)

  end subroutine s_parmatch_bld_wnxt

  function s_parmatch_aggregator_fmt() result(val)
    implicit none
    character(len=32)  :: val

    val = "Parallel Matching aggregation"
  end function s_parmatch_aggregator_fmt

  function amg_s_parmatch_aggregator_xt_desc() result(val)
    implicit none
    logical  :: val

    val = .true.
  end function amg_s_parmatch_aggregator_xt_desc

  function s_parmatch_aggregator_sizeof(ag) result(val)
    use psb_realloc_mod
    implicit none
    class(amg_s_parmatch_aggregator_type), intent(in)  :: ag
    integer(psb_epk_)  :: val

    val = 4
    val = val + psb_size(ag%w) + psb_size(ag%w_nxt)
    if (allocated(ag%ac))      val = val + ag%ac%sizeof()
    if (allocated(ag%base_a))  val = val + ag%base_a%sizeof()
    if (allocated(ag%prol))    val = val + ag%prol%sizeof()
    if (allocated(ag%restr))   val = val + ag%restr%sizeof()
    if (allocated(ag%desc_ac)) val = val + ag%desc_ac%sizeof()
    if (allocated(ag%base_desc)) val = val + ag%base_desc%sizeof()
    if (allocated(ag%desc_ax)) val = val + ag%desc_ax%sizeof()

  end function s_parmatch_aggregator_sizeof

  subroutine  s_parmatch_aggregator_descr(ag,parms,iout,info)
    implicit none
    class(amg_s_parmatch_aggregator_type), intent(in) :: ag
    type(amg_sml_parms), intent(in)   :: parms
    integer(psb_ipk_), intent(in)  :: iout
    integer(psb_ipk_), intent(out) :: info

    write(iout,*) 'Parallel Matching Aggregator'
    write(iout,*) '   Number of matching  sweeps: ',ag%n_sweeps
    write(iout,*) '   Matching algorithm         : MatchBoxP (PREIS)'
    write(iout,*) 'Aggregator object type: ',ag%fmt()
    call parms%mldescr(iout,info)

    return
  end subroutine s_parmatch_aggregator_descr

  function is_legal_malg(alg) result(val)
    logical :: val
    integer(psb_ipk_) :: alg

    val = (0==alg)
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


  subroutine  s_parmatch_aggregator_update_next(ag,agnext,info)
    use psb_realloc_mod
    implicit none
    class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
    class(amg_s_base_aggregator_type), target, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    !
    !
    select type(agnext)
    class is (amg_s_parmatch_aggregator_type)
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
      if (ag%unsmoothed_hierarchy) then
        agnext%unsmoothed_hierarchy = .true.
        call move_alloc(ag%rwdesc,agnext%base_desc)
        call move_alloc(ag%rwa,agnext%base_a)
      end if

    class default
      ! What should we do here?
    end select
    info = 0
  end subroutine s_parmatch_aggregator_update_next

  subroutine s_parmatch_aggr_csetc(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_s_parmatch_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    character(len=*), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='s_parmatch_aggr_cseti'
    info = psb_success_

    ! For now we ignore IDX

    select case(psb_toupper(trim(what)))
    case('PRMC_REPRODUCIBLE_MATCHING')
      select case(psb_toupper(trim(val)))
      case('F','FALSE')
        ag%reproducible_matching = .false.
      case('REPRODUCIBLE','TRUE','T')
        ag%reproducible_matching =.true.
      end select
    case('PRMC_NEED_SYMMETRIZE')
      select case(psb_toupper(trim(val)))
      case('FALSE','F')
        ag%need_symmetrize = .false.
      case('SYMMETRIZE','TRUE','T')
        ag%need_symmetrize =.true.
      end select
    case('PRMC_UNSMOOTHED_HIERARCHY')
      select case(psb_toupper(trim(val)))
      case('F','FALSE')
        ag%unsmoothed_hierarchy = .false.
      case('T','TRUE')
        ag%unsmoothed_hierarchy =.true.
      end select
    case default
      ! Do nothing
    end select
    return
  end subroutine s_parmatch_aggr_csetc

  subroutine s_parmatch_aggr_cseti(ag,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_s_parmatch_aggregator_type), intent(inout) :: ag
    character(len=*), intent(in)                  :: what
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act, iwhat
    character(len=20)  :: name='s_parmatch_aggr_cseti'
    info = psb_success_

    ! For now we ignore IDX

    select case(psb_toupper(trim(what)))
    case('PRMC_MATCH_ALG')
      ag%matching_alg=val
    case('PRMC_SWEEPS')
      ag%n_sweeps=val
    case('AGGR_SIZE')
      ag%orig_aggr_size = val
      ag%n_sweeps=max(1,ceiling(log(val*1.0)/log(2.0)))
    case('PRMC_MAX_CSIZE')
      ag%max_csize=val
    case('PRMC_MAX_NLEVELS')
      ag%max_nlevels=val
    case('PRMC_W_SIZE')
      call ag%bld_default_w(val)
    case('PRMC_REPRODUCIBLE_MATCHING')
      ag%reproducible_matching = (val == 1)
    case('PRMC_NEED_SYMMETRIZE')
      ag%need_symmetrize = (val == 1)
    case('PRMC_UNSMOOTHED_HIERARCHY')
      ag%unsmoothed_hierarchy = (val == 1)
    case default
      ! Do nothing
    end select
    return
  end subroutine s_parmatch_aggr_cseti

  subroutine s_parmatch_aggr_set_default(ag)

    Implicit None

    ! Arguments
    class(amg_s_parmatch_aggregator_type), intent(inout) :: ag
    character(len=20)  :: name='s_parmatch_aggr_set_default'
    call ag%amg_s_base_aggregator_type%default()
    ag%matching_alg  = 0
    ag%n_sweeps      = 1
    ag%jacobi_sweeps = 0
    ag%max_nlevels   = 36
    ag%max_csize     = -1
    !
    ! Apparently BootCMatch works better
    ! by keeping all entries
    !
    ag%do_clean_zeros = .false.

    return

  end subroutine s_parmatch_aggr_set_default

  subroutine  s_parmatch_aggregator_free(ag,info)
    use iso_c_binding
    implicit none
    class(amg_s_parmatch_aggregator_type), intent(inout) :: ag
    integer(psb_ipk_), intent(out)       :: info

    info = 0
    if ((info == 0).and.allocated(ag%w)) deallocate(ag%w,stat=info)
    if ((info == 0).and.allocated(ag%w_nxt)) deallocate(ag%w_nxt,stat=info)
    if ((info == 0).and.allocated(ag%prol)) then
      call ag%prol%free(); deallocate(ag%prol,stat=info)
    end if
    if ((info == 0).and.allocated(ag%restr)) then
      call ag%restr%free(); deallocate(ag%restr,stat=info)
    end if
    if ((info == 0).and.allocated(ag%ac)) then
      call ag%ac%free(); deallocate(ag%ac,stat=info)
    end if
    if ((info == 0).and.allocated(ag%base_a)) then
      call ag%base_a%free(); deallocate(ag%base_a,stat=info)
    end if
    if ((info == 0).and.allocated(ag%rwa)) then
      call ag%rwa%free(); deallocate(ag%rwa,stat=info)
    end if
    if ((info == 0).and.allocated(ag%desc_ac)) then
      call ag%desc_ac%free(info); deallocate(ag%desc_ac,stat=info)
    end if
    if ((info == 0).and.allocated(ag%desc_ax)) then
      call ag%desc_ax%free(info); deallocate(ag%desc_ax,stat=info)
    end if
    if ((info == 0).and.allocated(ag%base_desc)) then
      call ag%base_desc%free(info); deallocate(ag%base_desc,stat=info)
    end if
    if ((info == 0).and.allocated(ag%rwdesc)) then
      call ag%rwdesc%free(info); deallocate(ag%rwdesc,stat=info)
    end if

  end subroutine s_parmatch_aggregator_free

  subroutine  s_parmatch_aggregator_clone(ag,agnext,info)
    implicit none
    class(amg_s_parmatch_aggregator_type), intent(inout) :: ag
    class(amg_s_base_aggregator_type), allocatable, intent(inout) :: agnext
    integer(psb_ipk_), intent(out)       :: info

    info = 0
    if (allocated(agnext)) then
      call agnext%free(info)
      if (info == 0) deallocate(agnext,stat=info)
    end if
    if (info /= 0) return
    allocate(agnext,source=ag,stat=info)
    select type(agnext)
    class is (amg_s_parmatch_aggregator_type)
      call agnext%set_c_default_w()
    class default
      ! Should never ever get here
      info = -1
    end select
  end subroutine s_parmatch_aggregator_clone

  subroutine  amg_s_parmatch_aggregator_bld_map(ag,desc_a,desc_ac,ilaggr,nlaggr,&
       & op_restr,op_prol,map,info)
    use psb_base_mod
    implicit none
    class(amg_s_parmatch_aggregator_type), target, intent(inout) :: ag
    type(psb_desc_type), intent(in), target :: desc_a, desc_ac
    integer(psb_lpk_), intent(inout)        :: ilaggr(:), nlaggr(:)
    type(psb_sspmat_type), intent(inout)   :: op_prol, op_restr
    type(psb_slinmap_type), intent(out)     :: map
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='s_parmatch_aggregator_bld_map'

    call psb_erractionsave(err_act)
    !
    ! Copy the prolongation/restriction matrices into the descriptor map.
    !  op_restr => PR^T   i.e. restriction  operator
    !  op_prol => PR     i.e. prolongation operator
    !
    !  For parmatch have an explicit copy of the descriptors
    !
    if (allocated(ag%desc_ax)) then
!!$      write(0,*) 'Building linmap with ag%desc_ax ',ag%desc_ax%get_local_rows(),ag%desc_ax%get_local_cols(),&
!!$           & desc_ac%get_local_rows(),desc_ac%get_local_cols()
      map = psb_linmap(psb_map_gen_linear_,ag%desc_ax,&
           & desc_ac,op_restr,op_prol,ilaggr,nlaggr)
    else
      map = psb_linmap(psb_map_gen_linear_,desc_a,&
           & desc_ac,op_restr,op_prol,ilaggr,nlaggr)
    end if
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='sp_Free')
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine amg_s_parmatch_aggregator_bld_map
#endif
end module amg_s_parmatch_aggregator_mod

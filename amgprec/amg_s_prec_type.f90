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
! File: amg_s_prec_type.f90
!
! Module: amg_s_prec_type
!
!  This module defines:
!  - the amg_s_prec_type data structure containing the preconditioner and related
!    data structures;
!
!  It contains routines for
!  - Building and applying;
!  - checking if the preconditioner is correctly defined;
!  - printing a	description of the preconditioner;
!  - deallocating the preconditioner data structure.
!

module amg_s_prec_type

  use amg_base_prec_type
  use amg_s_base_solver_mod
  use amg_s_base_smoother_mod
  use amg_s_base_aggregator_mod
  use amg_s_onelev_mod
  use psb_base_mod, only : psb_erractionsave, psb_erractionrestore, &
       & psb_errstatus_fatal, psb_ctxt_type
  use psb_prec_mod, only : psb_sprec_type

  !
  ! Type: amg_sprec_type.
  !
  !  This is the data type containing all the information about the multilevel
  !  preconditioner ('d', 's', 'c' and 'z', according to the real/complex,
  !  single/double precision version of  AMG4PSBLAS).
  !  It consists of an array of 'one-level' intermediate data structures
  !  of type amg_sonelev_type, each containing the information needed to apply
  !  the smoothing and the coarse-space correction at a generic level. RT is the
  !  real data type, i.e. S for both S and C, and D for both D and Z.
  !
  !  type amg_sprec_type
  !    type(amg_sonelev_type), allocatable :: precv(:)
  !  end type amg_sprec_type
  !
  !  Note that the levels are numbered in increasing order starting from
  !  the level 1 as the finest one, and the number of levels is given by
  !  size(precv(:))  which is the id of the coarsest level.
  !  In the multigrid literature many authors number the levels in the opposite
  !  order, with level 0 being the id of the coarsest level.
  !
  !
  integer, parameter, private :: wv_size_=4

  type, extends(psb_sprec_type)        :: amg_sprec_type
    type(amg_saggr_data)                 :: ag_data
    !
    ! Number of outer sweeps. Sometimes  2 V-cycles may be better than 1 W-cycle.
    !
    integer(psb_ipk_)                  :: outer_sweeps = 1
    !
    ! Coarse solver requires some tricky checks, and for this we need to
    ! record the choice in the format given by the user,
    ! to keep track against what is put later in the multilevel array
    !
    integer(psb_ipk_)                  :: coarse_solver = -1

    !
    ! The multilevel hierarchy
    !
    type(amg_s_onelev_type), allocatable :: precv(:)
  contains
    procedure, pass(prec)               :: psb_s_apply2_vect => amg_s_apply2_vect
    procedure, pass(prec)               :: psb_s_apply1_vect => amg_s_apply1_vect
    procedure, pass(prec)               :: psb_s_apply2v => amg_s_apply2v
    procedure, pass(prec)               :: psb_s_apply1v => amg_s_apply1v
    procedure, pass(prec)               :: dump           => amg_s_dump
    procedure, pass(prec)               :: cnv            => amg_s_cnv
    procedure, pass(prec)               :: clone          => amg_s_clone
    procedure, pass(prec)               :: free           => amg_s_prec_free
    procedure, pass(prec)               :: allocate_wrk   => amg_s_allocate_wrk
    procedure, pass(prec)               :: free_wrk       => amg_s_free_wrk
    procedure, pass(prec)               :: is_allocated_wrk => amg_s_is_allocated_wrk
    procedure, pass(prec)               :: get_complexity => amg_s_get_compl
    procedure, pass(prec)               :: cmp_complexity => amg_s_cmp_compl
    procedure, pass(prec)               :: get_avg_cr => amg_s_get_avg_cr
    procedure, pass(prec)               :: cmp_avg_cr => amg_s_cmp_avg_cr
    procedure, pass(prec)               :: get_nlevs  => amg_s_get_nlevs
    procedure, pass(prec)               :: get_nzeros => amg_s_get_nzeros
    procedure, pass(prec)               :: sizeof => amg_sprec_sizeof
    procedure, pass(prec)               :: setsm  => amg_sprecsetsm
    procedure, pass(prec)               :: setsv  => amg_sprecsetsv
    procedure, pass(prec)               :: setag  => amg_sprecsetag
    procedure, pass(prec)               :: cseti  => amg_scprecseti
    procedure, pass(prec)               :: csetc  => amg_scprecsetc
    procedure, pass(prec)               :: csetr  => amg_scprecsetr
    generic, public                     :: set => setsm, setsv, setag 
    procedure, pass(prec)               :: get_smoother => amg_s_get_smootherp
    procedure, pass(prec)               :: get_solver   => amg_s_get_solverp
    procedure, pass(prec)               :: move_alloc   => s_prec_move_alloc
    procedure, pass(prec)               :: init         => amg_sprecinit
    procedure, pass(prec)               :: build        => amg_sprecbld
    procedure, pass(prec)               :: hierarchy_build   => amg_s_hierarchy_bld
    procedure, pass(prec)               :: hierarchy_rebuild => amg_s_hierarchy_rebld
    procedure, pass(prec)               :: smoothers_build   => amg_s_smoothers_bld
    procedure, pass(prec)               :: descr        =>  amg_sfile_prec_descr
  end type amg_sprec_type

  private :: amg_s_dump, amg_s_get_compl,  amg_s_cmp_compl,&
       & amg_s_get_avg_cr,  amg_s_cmp_avg_cr,&
       & amg_s_get_nzeros, amg_s_get_nlevs, s_prec_move_alloc


  !
  ! Interfaces to routines for checking the definition of the preconditioner,
  ! for printing its description and for deallocating its data structure
  !

  interface amg_precfree
    module procedure amg_sprecfree
  end interface


  interface amg_precdescr
    subroutine amg_sfile_prec_descr(prec,info,iout,root,verbosity)
      import :: amg_sprec_type, psb_ipk_
      implicit none
      ! Arguments
      class(amg_sprec_type), intent(in)      :: prec
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), intent(in), optional :: iout
      integer(psb_ipk_), intent(in), optional :: root
      integer(psb_ipk_), intent(in), optional :: verbosity
    end subroutine amg_sfile_prec_descr
  end interface

  interface amg_sizeof
    module procedure amg_sprec_sizeof
  end interface

  interface amg_precapply
    subroutine amg_sprecaply2_vect(prec,x,y,desc_data,info,trans,work)
      import :: psb_sspmat_type, psb_desc_type, &
           & psb_spk_, psb_s_vect_type, amg_sprec_type, psb_ipk_
      type(psb_desc_type),intent(in)      :: desc_data
      type(amg_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type),intent(inout) :: x
      type(psb_s_vect_type),intent(inout) :: y
      integer(psb_ipk_), intent(out)                :: info
      character(len=1), optional          :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine amg_sprecaply2_vect
    subroutine amg_sprecaply1_vect(prec,x,desc_data,info,trans,work)
      import :: psb_sspmat_type, psb_desc_type, &
           & psb_spk_, psb_s_vect_type, amg_sprec_type, psb_ipk_
      type(psb_desc_type),intent(in)      :: desc_data
      type(amg_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type),intent(inout) :: x
      integer(psb_ipk_), intent(out)                :: info
      character(len=1), optional          :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine amg_sprecaply1_vect
    subroutine amg_sprecaply(prec,x,y,desc_data,info,trans,work)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, amg_sprec_type, psb_ipk_
      type(psb_desc_type),intent(in)   :: desc_data
      type(amg_sprec_type), intent(inout) :: prec
      real(psb_spk_),intent(inout)     :: x(:)
      real(psb_spk_),intent(inout)     :: y(:)
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), optional       :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine amg_sprecaply
    subroutine amg_sprecaply1(prec,x,desc_data,info,trans)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, amg_sprec_type, psb_ipk_
      type(psb_desc_type),intent(in)   :: desc_data
      type(amg_sprec_type), intent(inout) :: prec
      real(psb_spk_),intent(inout)     :: x(:)
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), optional       :: trans
    end subroutine amg_sprecaply1
  end interface

  interface
    subroutine amg_sprecsetsm(prec,val,info,ilev,ilmax,pos)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, amg_s_base_smoother_type, psb_ipk_
      class(amg_sprec_type), target, intent(inout):: prec
      class(amg_s_base_smoother_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
      character(len=*), optional, intent(in)      :: pos
    end subroutine amg_sprecsetsm
    subroutine amg_sprecsetsv(prec,val,info,ilev,ilmax,pos)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, amg_s_base_solver_type, psb_ipk_
      class(amg_sprec_type), intent(inout)      :: prec
      class(amg_s_base_solver_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
      character(len=*), optional, intent(in)      :: pos
    end subroutine amg_sprecsetsv
    subroutine amg_sprecsetag(prec,val,info,ilev,ilmax,pos)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, amg_s_base_aggregator_type, psb_ipk_
      class(amg_sprec_type), intent(inout)      :: prec
      class(amg_s_base_aggregator_type), intent(in) :: val
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: ilev,ilmax
      character(len=*), optional, intent(in)      :: pos
    end subroutine amg_sprecsetag
    subroutine amg_scprecseti(prec,what,val,info,ilev,ilmax,pos,idx)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, psb_ipk_
      class(amg_sprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what
      integer(psb_ipk_), intent(in)            :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
      character(len=*), optional, intent(in)      :: pos
    end subroutine amg_scprecseti
    subroutine amg_scprecsetr(prec,what,val,info,ilev,ilmax,pos,idx)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, psb_ipk_
      class(amg_sprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what
      real(psb_spk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
      character(len=*), optional, intent(in)      :: pos
    end subroutine amg_scprecsetr
    subroutine amg_scprecsetc(prec,what,string,info,ilev,ilmax,pos,idx)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, psb_ipk_
      class(amg_sprec_type), intent(inout)   :: prec
      character(len=*), intent(in)             :: what
      character(len=*), intent(in)             :: string
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
      character(len=*), optional, intent(in)      :: pos
    end subroutine amg_scprecsetc
  end interface

  interface amg_precinit
    subroutine amg_sprecinit(ctxt,prec,ptype,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, psb_ipk_, psb_ctxt_type
      type(psb_ctxt_type), intent(in)        :: ctxt
      class(amg_sprec_type), intent(inout) :: prec
      character(len=*), intent(in)           :: ptype
      integer(psb_ipk_), intent(out)         :: info
    end subroutine amg_sprecinit
  end interface amg_precinit

  interface amg_precbld
    subroutine amg_sprecbld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & psb_s_base_sparse_mat, psb_s_base_vect_type, &
           & psb_i_base_vect_type, amg_sprec_type, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(amg_sprec_type), intent(inout), target       :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: amold
      class(psb_s_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      !      character, intent(in),optional             :: upd
    end subroutine amg_sprecbld
  end interface amg_precbld

  interface amg_hierarchy_bld
    subroutine amg_s_hierarchy_bld(a,desc_a,prec,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(amg_sprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      !      character, intent(in),optional             :: upd
    end subroutine amg_s_hierarchy_bld
  end interface amg_hierarchy_bld

  interface amg_hierarchy_rebld
    subroutine amg_s_hierarchy_rebld(a,desc_a,prec,info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & amg_sprec_type, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(amg_sprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      !      character, intent(in),optional             :: upd
    end subroutine amg_s_hierarchy_rebld
  end interface amg_hierarchy_rebld

  interface amg_smoothers_bld
    subroutine amg_s_smoothers_bld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, &
           & psb_s_base_sparse_mat, psb_s_base_vect_type, &
           & psb_i_base_vect_type, amg_sprec_type, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(in), target          :: a
      type(psb_desc_type), intent(inout), target           :: desc_a
      class(amg_sprec_type), intent(inout), target        :: prec
      integer(psb_ipk_), intent(out)                       :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: amold
      class(psb_s_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
      !      character, intent(in),optional             :: upd
    end subroutine amg_s_smoothers_bld
  end interface amg_smoothers_bld

contains
  !
  ! Function returning a pointer to the smoother
  !
  function amg_s_get_smootherp(prec,ilev) result(val)
    implicit none
    class(amg_sprec_type), target, intent(in) :: prec
    integer(psb_ipk_), optional                 :: ilev
    class(amg_s_base_smoother_type), pointer  :: val
    integer(psb_ipk_)        :: ilev_

    val => null()
    if (present(ilev)) then
      ilev_ = ilev
    else
      ! What is a good default?
      ilev_ = 1
    end if
    if (allocated(prec%precv)) then
      if ((1<=ilev_).and.(ilev_<=size(prec%precv))) then
        if (allocated(prec%precv(ilev_)%sm)) then
          val => prec%precv(ilev_)%sm
        end if
      end if
    end if
  end function amg_s_get_smootherp
  !
  ! Function returning a pointer to the solver
  !
  function amg_s_get_solverp(prec,ilev) result(val)
    implicit none
    class(amg_sprec_type), target, intent(in) :: prec
    integer(psb_ipk_), optional                 :: ilev
    class(amg_s_base_solver_type), pointer  :: val
    integer(psb_ipk_)        :: ilev_

    val => null()
    if (present(ilev)) then
      ilev_ = ilev
    else
      ! What is a good default?
      ilev_ = 1
    end if
    if (allocated(prec%precv)) then
      if ((1<=ilev_).and.(ilev_<=size(prec%precv))) then
        if (allocated(prec%precv(ilev_)%sm)) then
          if (allocated(prec%precv(ilev_)%sm%sv)) then
            val => prec%precv(ilev_)%sm%sv
          end if
        end if
      end if
    end if
  end function amg_s_get_solverp
  !
  ! Function returning the size of the precv(:) array
  !
  function amg_s_get_nlevs(prec) result(val)
    implicit none
    class(amg_sprec_type), intent(in) :: prec
    integer(psb_ipk_) :: val
    val = 0
    if (allocated(prec%precv)) then
      val = size(prec%precv)
    end if
  end function amg_s_get_nlevs
  !
  ! Function returning the size of the amg_prec_type data structure
  ! in bytes or in number of nonzeros of the operator(s) involved.
  !
  function amg_s_get_nzeros(prec) result(val)
    implicit none
    class(amg_sprec_type), intent(in) :: prec
    integer(psb_epk_) :: val
    integer(psb_ipk_)        :: i
    val = 0
    if (allocated(prec%precv)) then
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%get_nzeros()
      end do
    end if
  end function amg_s_get_nzeros

  function amg_sprec_sizeof(prec, global) result(val)
    implicit none
    class(amg_sprec_type), intent(in) :: prec
    logical, intent(in), optional :: global
    integer(psb_epk_) :: val    
    integer(psb_ipk_)        :: i
    type(psb_ctxt_type) :: ctxt
    
    logical :: global_

    if (present(global)) then
      global_ = global
    else
      global_ = .false.
    end if
    
    val = 0
    val = val + psb_sizeof_ip
    if (allocated(prec%precv)) then
      do i=1, size(prec%precv)
        val = val + prec%precv(i)%sizeof()
      end do
    end if
    if (global_) then
      ctxt = prec%ctxt
      call psb_sum(ctxt,val)
    end if

  end function amg_sprec_sizeof

  !
  ! Operator complexity: ratio of total number
  ! of nonzeros in the aggregated matrices at the
  ! various level to the nonzeroes at the fine level
  ! (original matrix)
  !

  function amg_s_get_compl(prec) result(val)
    implicit none
    class(amg_sprec_type), intent(in) :: prec
    real(psb_spk_)  :: val

    val = prec%ag_data%op_complexity

  end function amg_s_get_compl

  subroutine amg_s_cmp_compl(prec)

    implicit none
    class(amg_sprec_type), intent(inout) :: prec

    real(psb_spk_) :: num, den, nmin
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: il

    num = -sone
    den = sone
    ctxt = prec%ctxt
    if (allocated(prec%precv)) then
      il  = 1
      num = prec%precv(il)%base_a%get_nzeros()
      if (num >= szero) then
        den = num
        do il=2,size(prec%precv)
          num = num + max(0,prec%precv(il)%base_a%get_nzeros())
        end do
      end if
    end if
    nmin = num
    call psb_min(ctxt,nmin)
    if (nmin < szero) then
      num = szero
      den = sone
    else
      call psb_sum(ctxt,num)
      call psb_sum(ctxt,den)
    end if
    prec%ag_data%op_complexity = num/den
  end subroutine amg_s_cmp_compl

  !
  ! Average coarsening ratio
  !

  function amg_s_get_avg_cr(prec) result(val)
    implicit none
    class(amg_sprec_type), intent(in) :: prec
    real(psb_spk_)  :: val

    val = prec%ag_data%avg_cr

  end function amg_s_get_avg_cr

  subroutine amg_s_cmp_avg_cr(prec)

    implicit none
    class(amg_sprec_type), intent(inout) :: prec

    real(psb_spk_)    :: avgcr
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: il, nl, iam, np


    avgcr = szero
    ctxt = prec%ctxt
    call psb_info(ctxt,iam,np)
    if (allocated(prec%precv)) then
      nl = size(prec%precv)
      do il=2,nl
        avgcr = avgcr + max(szero,prec%precv(il)%szratio)
      end do
      avgcr = avgcr / (nl-1)
    end if
    call psb_sum(ctxt,avgcr)
    prec%ag_data%avg_cr = avgcr/np
  end subroutine amg_s_cmp_avg_cr

  !
  ! Subroutines: amg_Tprec_free
  ! Version: real
  !
  !  These routines deallocate the amg_Tprec_type data structures.
  !
  ! Arguments:
  !  p       -  type(amg_Tprec_type), input.
  !             The data structure to be deallocated.
  !  info    -  integer, output.
  !             error code.
  !
  subroutine amg_sprecfree(p,info)

    implicit none

    ! Arguments
    type(amg_sprec_type), intent(inout) :: p
    integer(psb_ipk_), intent(out)        :: info

    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i
    character(len=20)   :: name

    info=psb_success_
    name = 'amg_sprecfree'
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; return
    end if

    me=-1

    call p%free(info)


    return

  end subroutine amg_sprecfree

  subroutine amg_s_prec_free(prec,info)

    implicit none

    ! Arguments
    class(amg_sprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info

    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i
    character(len=20)   :: name

    info=psb_success_
    name = 'amg_sprecfree'
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if

    me=-1
    if (allocated(prec%precv)) then
      do i=1,size(prec%precv)
        call prec%precv(i)%free(info)
      end do
      deallocate(prec%precv,stat=info)
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_s_prec_free



  !
  ! Top level methods.
  !
  subroutine amg_s_apply2_vect(prec,x,y,desc_data,info,trans,work)
    implicit none
    type(psb_desc_type),intent(in)        :: desc_data
    class(amg_sprec_type), intent(inout)  :: prec
    type(psb_s_vect_type),intent(inout)   :: x
    type(psb_s_vect_type),intent(inout)   :: y
    integer(psb_ipk_), intent(out)          :: info
    character(len=1), optional              :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec)
    type is (amg_sprec_type)
      call amg_precapply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_s_apply2_vect

  subroutine amg_s_apply1_vect(prec,x,desc_data,info,trans,work)
    implicit none
    type(psb_desc_type),intent(in)          :: desc_data
    class(amg_sprec_type), intent(inout)  :: prec
    type(psb_s_vect_type),intent(inout)   :: x
    integer(psb_ipk_), intent(out)          :: info
    character(len=1), optional            :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec)
    type is (amg_sprec_type)
      call amg_precapply(prec,x,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_s_apply1_vect


  subroutine amg_s_apply2v(prec,x,y,desc_data,info,trans,work)
    implicit none
    type(psb_desc_type),intent(in)    :: desc_data
    class(amg_sprec_type), intent(inout) :: prec
    real(psb_spk_),intent(inout)      :: x(:)
    real(psb_spk_),intent(inout)      :: y(:)
    integer(psb_ipk_), intent(out)     :: info
    character(len=1), optional        :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec)
    type is (amg_sprec_type)
      call amg_precapply(prec,x,y,desc_data,info,trans,work)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_s_apply2v

  subroutine amg_s_apply1v(prec,x,desc_data,info,trans)
    implicit none
    type(psb_desc_type),intent(in)    :: desc_data
    class(amg_sprec_type), intent(inout) :: prec
    real(psb_spk_),intent(inout)      :: x(:)
    integer(psb_ipk_), intent(out)     :: info
    character(len=1), optional         :: trans
    Integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_prec_apply'

    call psb_erractionsave(err_act)

    select type(prec)
    type is (amg_sprec_type)
      call amg_precapply(prec,x,desc_data,info,trans)
    class default
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
  return

  end subroutine amg_s_apply1v


  subroutine amg_s_dump(prec,info,istart,iend,iproc,prefix,head,&
       & ac,rp,smoother,solver,tprol,&
       & global_num)

    implicit none
    class(amg_sprec_type), intent(in)     :: prec
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_), intent(in), optional :: istart, iend, iproc
    character(len=*), intent(in), optional  :: prefix, head
    logical, optional, intent(in)    :: smoother, solver,ac, rp, tprol, global_num
    integer(psb_ipk_)   :: i, j, il1, iln, lev
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: iam, np, iproc_
    character(len=80)   :: prefix_
    character(len=120)  :: fname ! len should be at least 20 more than
    !  len of prefix_

    info = 0
    ctxt = prec%ctxt
    call psb_info(ctxt,iam,np)
    iln = size(prec%precv)
    if (present(istart)) then
      il1 = max(1,istart)
    else
      il1 = min(2,iln)
    end if
    if (present(iend)) then
      iln = min(iln, iend)
    end if
    iproc_ = -1
    if (present(iproc)) then
      iproc_ = iproc
    end if

    if ((iproc_ == -1).or.(iproc_==iam)) then
      do lev=il1, iln
        call prec%precv(lev)%dump(lev,info,prefix=prefix,head=head,&
             & ac=ac,smoother=smoother,solver=solver,rp=rp,tprol=tprol, &
             & global_num=global_num)
      end do
    end if
  end subroutine amg_s_dump

  subroutine amg_s_cnv(prec,info,amold,vmold,imold)

    implicit none
    class(amg_sprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)       :: info
    class(psb_s_base_sparse_mat), intent(in), optional :: amold
    class(psb_s_base_vect_type), intent(in), optional  :: vmold
    class(psb_i_base_vect_type), intent(in), optional  :: imold

    integer(psb_ipk_) :: i

    info = psb_success_
    if (allocated(prec%precv)) then
      do i=1,size(prec%precv)
        if (info == psb_success_ ) &
             & call prec%precv(i)%cnv(info,amold=amold,vmold=vmold,imold=imold)
      end do
    end if

  end subroutine amg_s_cnv

  subroutine amg_s_clone(prec,precout,info)

    implicit none
    class(amg_sprec_type), intent(inout) :: prec
    class(psb_sprec_type), intent(inout) :: precout
    integer(psb_ipk_), intent(out)       :: info

    call precout%free(info)
    if (info == 0) call amg_s_inner_clone(prec,precout,info)

  end subroutine amg_s_clone

  subroutine amg_s_inner_clone(prec,precout,info)

    implicit none
    class(amg_sprec_type), intent(inout)         :: prec
    class(psb_sprec_type), target, intent(inout) :: precout
    integer(psb_ipk_), intent(out)             :: info
    ! Local vars
    integer(psb_ipk_)   :: i, j, ln, lev
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: iam, np

    info = psb_success_
    select type(pout => precout)
    class is (amg_sprec_type)
      pout%ctxt          = prec%ctxt
      pout%ag_data       = prec%ag_data
      pout%outer_sweeps  = prec%outer_sweeps
      if (allocated(prec%precv)) then
        ln = size(prec%precv)
        allocate(pout%precv(ln),stat=info)
        if (info /= psb_success_) goto 9999
        if (ln >= 1) then
          call prec%precv(1)%clone(pout%precv(1),info)
        end if
        do lev=2, ln
          if (info /= psb_success_) exit
          call prec%precv(lev)%clone(pout%precv(lev),info)
          if (info == psb_success_) then
            pout%precv(lev)%base_a       => pout%precv(lev)%ac
            pout%precv(lev)%base_desc    => pout%precv(lev)%desc_ac
            pout%precv(lev)%linmap%p_desc_U => pout%precv(lev-1)%base_desc
            pout%precv(lev)%linmap%p_desc_V => pout%precv(lev)%base_desc
          end if
        end do
      end if
      if (allocated(prec%precv(1)%wrk)) &
           & call pout%allocate_wrk(info,vmold=prec%precv(1)%wrk%vx2l%v)

    class default
      write(0,*) 'Error: wrong out type'
      info = psb_err_invalid_input_
    end select
9999 continue
  end subroutine amg_s_inner_clone

  subroutine s_prec_move_alloc(prec, b,info)
    use psb_base_mod
    implicit none
    class(amg_sprec_type), intent(inout) :: prec
    class(amg_sprec_type), intent(inout), target :: b
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i

    if (same_type_as(prec,b)) then
      if (allocated(b%precv)) then
        ! This might not be required if FINAL procedures are available.
        call b%free(info)
        if (info /= psb_success_) then
          !?????
!!$        return
        endif
      end if
      b%ctxt         = prec%ctxt
      b%ag_data       = prec%ag_data
      b%outer_sweeps  = prec%outer_sweeps

      call move_alloc(prec%precv,b%precv)
      ! Fix the pointers except on level 1.
      do i=2, size(b%precv)
        b%precv(i)%base_a    => b%precv(i)%ac
        b%precv(i)%base_desc => b%precv(i)%desc_ac
        b%precv(i)%linmap%p_desc_U => b%precv(i-1)%base_desc
        b%precv(i)%linmap%p_desc_V => b%precv(i)%base_desc
      end do

    else
      write(0,*) 'Warning: PREC%move_alloc onto different type?'
      info = psb_err_internal_error_
    end if
  end subroutine s_prec_move_alloc

  subroutine amg_s_allocate_wrk(prec,info,vmold,desc)
    use psb_base_mod
    implicit none

    ! Arguments
    class(amg_sprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info
    class(psb_s_base_vect_type), intent(in), optional  :: vmold
    !
    ! In MLD the DESC optional argument is ignored, since
    ! the necessary info is contained in the various entries of the
    ! PRECV component.
    type(psb_desc_type), intent(in), optional  :: desc

    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i,j,level,nlev, nc2l
    character(len=20)   :: name

    info=psb_success_
    name = 'amg_s_allocate_wrk'
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if
    nlev   = size(prec%precv)
    level = 1
    do level = 1, nlev
      call prec%precv(level)%allocate_wrk(info,vmold=vmold)
      if (psb_errstatus_fatal()) then
        nc2l = prec%precv(level)%base_desc%get_local_cols()
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/2*nc2l/), a_err='real(psb_spk_)')
        goto 9999
      end if
    end do

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_s_allocate_wrk

  subroutine amg_s_free_wrk(prec,info)
    use psb_base_mod
    implicit none

    ! Arguments
    class(amg_sprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info

    ! Local variables
    integer(psb_ipk_)   :: me,err_act,i,j,level, nlev, nc2l
    character(len=20)   :: name

    info=psb_success_
    name = 'amg_s_free_wrk'
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if

    if (allocated(prec%precv)) then
      nlev   = size(prec%precv)
      do level = 1, nlev
        call prec%precv(level)%free_wrk(info)
      end do
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_s_free_wrk

  function amg_s_is_allocated_wrk(prec) result(res)
    use psb_base_mod
    implicit none

    ! Arguments
    class(amg_sprec_type), intent(in) :: prec
    logical :: res

    res = .false.
    if (.not.allocated(prec%precv)) return
    res = allocated(prec%precv(1)%wrk)

  end function amg_s_is_allocated_wrk

end module amg_s_prec_type

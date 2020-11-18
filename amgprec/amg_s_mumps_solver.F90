
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
!  Current version of this file contributed by:
!        Ambra Abdullahi Hassan 
!
!
! File: amg_s_mumps_solver_mod.f90
!
! Module: amg_s_mumps_solver_mod
!
!  This module defines: 
!  - the amg_s_mumps_solver_type data structure containing the ingredients
!    to interface with the MUMPS package. 
!    1. The factorization can be either restricted to the diagonal block of the
!       current image or distributed (and thus exact). 
!
module amg_s_mumps_solver
  use amg_s_base_solver_mod
#if defined(HAVE_MUMPS_) && defined(HAVE_MUMPS_MODULES_)
  use smumps_struc_def
#endif
#if defined(HAVE_MUMPS_) && defined(HAVE_MUMPS_INCLUDES_)
  include 'smumps_struc.h'
#endif  
  
  
  type :: amg_s_mumps_icntl_item
    integer(psb_ipk_), allocatable :: item
  end type amg_s_mumps_icntl_item
  type :: amg_s_mumps_rcntl_item
    real(psb_spk_), allocatable :: item
  end type amg_s_mumps_rcntl_item

  type, extends(amg_s_base_solver_type) :: amg_s_mumps_solver_type
#if defined(HAVE_MUMPS_)
    type(smumps_struc), allocatable  :: id
#else
    integer, allocatable :: id
#endif
    type(amg_s_mumps_icntl_item), allocatable :: icntl(:)
    type(amg_s_mumps_rcntl_item), allocatable :: rcntl(:)
    !
    ! Controls to be set before MUMPS instantiation:
    !
    ! IPAR(1) : MUMPS_LOC_GLOB   0==amg_local_solver_: LOCAL   1==amg_global_solver_: GLOBAL
    ! IPAR(2) : MUMPS_PRINT_ERR  print verbosity (see MUMPS)
    ! IPAR(3) : MUMPS_SYM        0: non-symmetric   2: symmetric
    integer(psb_ipk_), dimension(3) :: ipar
      type(psb_ctxt_type), allocatable  :: local_ctxt
    logical                         :: built = .false.
  contains
    procedure, pass(sv) :: build   => s_mumps_solver_bld
    procedure, pass(sv) :: apply_a => s_mumps_solver_apply
    procedure, pass(sv) :: apply_v => s_mumps_solver_apply_vect
    procedure, pass(sv) :: clone_settings  => s_mumps_solver_clone_settings
    procedure, pass(sv) :: clear_data  => s_mumps_solver_clear_data
    procedure, pass(sv) :: free    => s_mumps_solver_free
    procedure, pass(sv) :: descr   => s_mumps_solver_descr
    procedure, pass(sv) :: sizeof  => s_mumps_solver_sizeof
    procedure, pass(sv) :: csetc   => s_mumps_solver_csetc
    procedure, pass(sv) :: cseti   => s_mumps_solver_cseti
    procedure, pass(sv) :: csetr   => s_mumps_solver_csetr
    procedure, pass(sv) :: default => s_mumps_solver_default
    procedure, nopass   :: get_fmt => s_mumps_solver_get_fmt
    procedure, nopass   :: get_id  => s_mumps_solver_get_id
    procedure, pass(sv) :: is_global => s_mumps_solver_is_global
    final               :: s_mumps_solver_finalize
  end type amg_s_mumps_solver_type


  private :: s_mumps_solver_bld, s_mumps_solver_apply, &
       &  s_mumps_solver_free,   s_mumps_solver_descr, &
       &  s_mumps_solver_sizeof, s_mumps_solver_apply_vect,&
       &  s_mumps_solver_cseti, s_mumps_solver_csetr,   &
       &  s_mumps_solver_csetc,  s_mumps_solver_clear_data,  &
       &  s_mumps_solver_default, s_mumps_solver_get_fmt, &
       &  s_mumps_solver_clone_settings, &
       &  s_mumps_solver_get_id, s_mumps_solver_is_global
  private :: s_mumps_solver_finalize

  interface 
    subroutine s_mumps_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, amg_s_mumps_solver_type, psb_s_vect_type, psb_dpk_, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_s_mumps_solver_type), intent(inout) :: sv
      type(psb_s_vect_type),intent(inout)  :: x
      type(psb_s_vect_type),intent(inout)  :: y
      real(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_spk_),target, intent(inout) :: work(:)
      type(psb_s_vect_type),intent(inout) :: wv(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_s_vect_type),intent(inout), optional   :: initu
    end subroutine s_mumps_solver_apply_vect
  end interface

  interface
    subroutine s_mumps_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info,init,initu)
      import :: psb_desc_type, amg_s_mumps_solver_type, psb_s_vect_type, psb_dpk_, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_s_mumps_solver_type), intent(inout) :: sv
      real(psb_spk_),intent(inout)         :: x(:)
      real(psb_spk_),intent(inout)         :: y(:)
      real(psb_spk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_spk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)                 :: info
      character, intent(in), optional       :: init
      real(psb_spk_),intent(inout), optional :: initu(:)
    end subroutine s_mumps_solver_apply
  end interface

  interface
    subroutine s_mumps_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)

      import :: psb_desc_type, amg_s_mumps_solver_type, psb_s_vect_type, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type

      Implicit None

      ! Arguments
      type(psb_sspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(amg_s_mumps_solver_type), intent(inout)       :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_sspmat_type), intent(in), target, optional :: b
      class(psb_s_base_sparse_mat), intent(in), optional  :: amold
      class(psb_s_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine s_mumps_solver_bld
  end interface

contains

  subroutine s_mumps_solver_clone_settings(sv,svout,info)

    use psb_base_mod
    Implicit None
    ! Arguments
    class(amg_s_mumps_solver_type), intent(inout) :: sv
    class(amg_s_base_solver_type), intent(inout) :: svout
    integer(psb_ipk_), intent(out)                 :: info
    integer(psb_ipk_) :: k,err_act
    character(len=20) :: name='s_mumps_solver_clone_settings'

    info = 0

#if defined(HAVE_MUMPS_)
    
    call psb_erractionsave(err_act)

    select type(svout)
    class is(amg_s_mumps_solver_type)
      svout%ipar(:) = sv%ipar(:)
      svout%built = .false.
      if (allocated(svout%icntl)) deallocate(svout%icntl,stat=info)
      if (info == 0) allocate(svout%icntl(amg_mumps_icntl_size),stat=info)
      if (info == 0) then
        do k=1,amg_mumps_icntl_size
          call psb_safe_ab_cpy(sv%icntl(k)%item,svout%icntl(k)%item,info)
        end do
      end if
      
      if (allocated(svout%rcntl)) deallocate(svout%rcntl,stat=info)
      if (info == 0) allocate(svout%rcntl(amg_mumps_rcntl_size),stat=info)
      if (info == 0) then
        do k=1,amg_mumps_rcntl_size
          call psb_safe_ab_cpy(sv%rcntl(k)%item,svout%rcntl(k)%item,info)
        end do
      end if

    class default    
      info = psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999 
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
#endif
  end subroutine s_mumps_solver_clone_settings

  subroutine s_mumps_solver_clear_data(sv,info)
    use psb_base_mod, only : psb_exit
    Implicit None

    ! Arguments
    class(amg_s_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    Integer(psb_ipk_) :: err_act
    character(len=20)  :: name='s_mumps_solver_clear_data'

    info = 0
#if defined(HAVE_MUMPS_)
    call psb_erractionsave(err_act)
    if (allocated(sv%id)) then      
      if (sv%built) then 
        sv%id%job = -2
        call smumps(sv%id)
        info = sv%id%infog(1)
        if (info /= psb_success_) goto 9999
      end if
      deallocate(sv%id, stat=info)
      if (allocated(sv%local_ctxt)) then
        call psb_exit(sv%local_ctxt,close=.false.)
        deallocate(sv%local_ctxt,stat=info)
      end if
      sv%built=.false.
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
#endif
  end subroutine s_mumps_solver_clear_data

  subroutine s_mumps_solver_free(sv,info)
    use psb_base_mod, only : psb_exit
    Implicit None

    ! Arguments
    class(amg_s_mumps_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    Integer(psb_ipk_) :: err_act
    character(len=20)  :: name='s_mumps_solver_free'

    info = 0
#if defined(HAVE_MUMPS_)
    call psb_erractionsave(err_act)
    call sv%clear_data(info)
    if ((info == 0).and.allocated(sv%icntl)) deallocate(sv%icntl,stat=info)
    if ((info == 0).and.allocated(sv%rcntl)) deallocate(sv%rcntl,stat=info)

    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
#endif
  end subroutine s_mumps_solver_free

subroutine s_mumps_solver_finalize(sv)

  Implicit None

  ! Arguments
  type(amg_s_mumps_solver_type), intent(inout) :: sv 
  integer(psb_ipk_) :: info
  Integer(psb_ipk_) :: err_act
  character(len=20) :: name='s_mumps_solver_finalize'

  call sv%free(info) 

  return

end subroutine s_mumps_solver_finalize

subroutine s_mumps_solver_descr(sv,info,iout,coarse)

  Implicit None

  ! Arguments
  class(amg_s_mumps_solver_type), intent(in) :: sv
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_), intent(in), optional            :: iout
  logical, intent(in), optional       :: coarse

  ! Local variables
  integer(psb_ipk_)   :: err_act
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  character(len=20), parameter :: name='amg_z_mumps_solver_descr'
  integer(psb_ipk_) :: iout_

  call psb_erractionsave(err_act)
  info = psb_success_
  if (present(iout)) then 
    iout_ = iout 
  else
    iout_ = psb_out_unit
  endif

  write(iout_,*) '  MUMPS  Solver. '

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine s_mumps_solver_descr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  WARNING: OTHER PARAMETERS OF MUMPS COULD BE ADDED.                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine s_mumps_solver_csetc(sv,what,val,info,idx)

  Implicit None

  ! Arguments
  class(amg_s_mumps_solver_type), intent(inout) :: sv
  character(len=*), intent(in)                  :: what
  character(len=*), intent(in)                  :: val
  integer(psb_ipk_), intent(out)                :: info
  integer(psb_ipk_), intent(in), optional       :: idx
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='s_mumps_solver_csetc'

  info = psb_success_
  call psb_erractionsave(err_act)


  select case(psb_toupper(trim(what)))
#if defined(HAVE_MUMPS_)
  case('MUMPS_LOC_GLOB')
    sv%ipar(1) = sv%stringval(psb_toupper(trim(val)))
#endif
  case default
    call sv%amg_s_base_solver_type%set(what,val,info,idx=idx)
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine s_mumps_solver_csetc


subroutine s_mumps_solver_cseti(sv,what,val,info,idx)

  Implicit None

  ! Arguments
  class(amg_s_mumps_solver_type), intent(inout) :: sv
  character(len=*), intent(in)                  :: what
  integer(psb_ipk_), intent(in)                 :: val
  integer(psb_ipk_), intent(out)                :: info
  integer(psb_ipk_), intent(in), optional       :: idx
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='s_mumps_solver_cseti'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(what))
#if defined(HAVE_MUMPS_)
  case('MUMPS_LOC_GLOB')
    sv%ipar(1) = val
  case('MUMPS_PRINT_ERR')
    sv%ipar(2) = val
  case('MUMPS_SYM')
    sv%ipar(3) = val 
  case('MUMPS_IPAR_ENTRY')
    if(present(idx)) then
      ! Note: this will allocate %item
      sv%icntl(idx)%item = val
    end if
#endif
  case default
    call sv%amg_s_base_solver_type%set(what,val,info,idx=idx)
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine s_mumps_solver_cseti

subroutine s_mumps_solver_csetr(sv,what,val,info,idx)

  Implicit None

  ! Arguments
  class(amg_s_mumps_solver_type), intent(inout) :: sv
  character(len=*), intent(in)                  :: what
  real(psb_spk_), intent(in)                 :: val
  integer(psb_ipk_), intent(out)                :: info
  integer(psb_ipk_), intent(in), optional       :: idx
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='s_mumps_solver_csetr'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(what))
#if defined(HAVE_MUMPS_)
  case('MUMPS_RPAR_ENTRY')
    if(present(idx)) then 
      ! Note: this will allocate %item
      sv%rcntl(idx)%item = val
    end if
#endif
  case default
    call sv%amg_s_base_solver_type%set(what,val,info,idx=idx)
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine s_mumps_solver_csetr

!!NOTE: BY DEFAULT BLR is activated with a dropping parameter to 1d-4       !!
subroutine s_mumps_solver_default(sv)

  Implicit none

  !Argument
  class(amg_s_mumps_solver_type),intent(inout) :: sv
  integer(psb_ipk_) :: info
  integer(psb_ipk_)  :: err_act,ictx,icomm
  character(len=20)  :: name='s_mumps_default'

  info = psb_success_
  call psb_erractionsave(err_act)

#if defined(HAVE_MUMPS_)
  if (.not.allocated(sv%id)) then 
    allocate(sv%id,stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name,a_err='amg_smumps_default')
      goto 9999
    end if
    sv%built=.false.
  end if
  if (.not.allocated(sv%icntl)) then
    allocate(sv%icntl(amg_mumps_icntl_size),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name,a_err='amg_smumps_default')
      goto 9999
    end if
  end if
  if (.not.allocated(sv%rcntl)) then
    allocate(sv%rcntl(amg_mumps_rcntl_size),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name,a_err='amg_smumps_default')
      goto 9999
    end if
  end if
  ! INSTANTIATION OF sv%id needed to set parmater but mpi communicator needed
  ! sv%id%job = -1
  ! sv%id%par=1
  ! call dmumps(sv%id)
  sv%ipar    = 0
  sv%ipar(1) = amg_global_solver_
  !sv%ipar(10)=6
  !sv%ipar(11)=0
  !sv%ipar(12)=6

#endif
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine s_mumps_solver_default

function s_mumps_solver_sizeof(sv) result(val)

  implicit none 
  ! Arguments
  class(amg_s_mumps_solver_type), intent(in) :: sv
  integer(psb_epk_) :: val
  integer             :: i
#if defined(HAVE_MUMPS_)
  val = (sv%id%INFOG(22)+sv%id%INFOG(32))*1d+6
#else
  val = 0 
#endif
  ! val = 2*psb_sizeof_ip + psb_sizeof_dp
  ! val = val + sv%symbsize
  ! val = val + sv%numsize
  return
end function s_mumps_solver_sizeof

function s_mumps_solver_get_fmt() result(val)
  implicit none 
  character(len=32)  :: val

  val = "MUMPS solver"
end function s_mumps_solver_get_fmt

function s_mumps_solver_get_id() result(val)
  implicit none 
  integer(psb_ipk_)  :: val

  val = amg_mumps_
end function s_mumps_solver_get_id


function s_mumps_solver_is_global(sv) result(val)
  implicit none 
  class(amg_s_mumps_solver_type), intent(in) :: sv
  logical  :: val

  val =  (sv%ipar(1) == amg_global_solver_ )
end function s_mumps_solver_is_global

end module amg_s_mumps_solver


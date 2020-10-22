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
!
! File: amg_d_ilu_solver_mod.f90
!
! Module: amg_d_ilu_solver_mod
!
!  This module defines: 
!  - the amg_d_ilu_solver_type data structure containing the ingredients
!    for a local Incomplete LU factorization.
!    1. The factorization is always restricted to the diagonal block of the
!       current image (coherently with the definition of a SOLVER as a local
!       object)
!    2. The code provides support for both pattern-based ILU(K) and
!       threshold base ILU(T,L)
!    3. The diagonal is stored separately, so strictly speaking this is
!       an incomplete LDU factorization; 
!    4. The application phase is shared among all variants;
!
!
module amg_d_ilu_solver

  use amg_base_prec_type, only : amg_fact_names
  use amg_d_base_solver_mod
  use psb_d_ilu_fact_mod

  type, extends(amg_d_base_solver_type) :: amg_d_ilu_solver_type
    type(psb_dspmat_type)      :: l, u
    real(psb_dpk_), allocatable :: d(:)
    type(psb_d_vect_type)      :: dv
    integer(psb_ipk_)            :: fact_type, fill_in
    real(psb_dpk_)                :: thresh
  contains
    procedure, pass(sv) :: dump    => amg_d_ilu_solver_dmp
    procedure, pass(sv) :: check   => d_ilu_solver_check
    procedure, pass(sv) :: clone   => amg_d_ilu_solver_clone
    procedure, pass(sv) :: clone_settings => amg_d_ilu_solver_clone_settings
    procedure, pass(sv) :: clear_data     => amg_d_ilu_solver_clear_data
    procedure, pass(sv) :: build   => amg_d_ilu_solver_bld
    procedure, pass(sv) :: cnv     => amg_d_ilu_solver_cnv
    procedure, pass(sv) :: apply_v => amg_d_ilu_solver_apply_vect
    procedure, pass(sv) :: apply_a => amg_d_ilu_solver_apply
    procedure, pass(sv) :: free    => d_ilu_solver_free
    procedure, pass(sv) :: cseti   => d_ilu_solver_cseti
    procedure, pass(sv) :: csetc   => d_ilu_solver_csetc
    procedure, pass(sv) :: csetr   => d_ilu_solver_csetr
    procedure, pass(sv) :: descr   => d_ilu_solver_descr
    procedure, pass(sv) :: default => d_ilu_solver_default
    procedure, pass(sv) :: sizeof  => d_ilu_solver_sizeof
    procedure, pass(sv) :: get_nzeros => d_ilu_solver_get_nzeros
    procedure, nopass   :: get_wrksz => d_ilu_solver_get_wrksize
    procedure, nopass   :: get_fmt    => d_ilu_solver_get_fmt
    procedure, nopass   :: get_id     => d_ilu_solver_get_id
  end type amg_d_ilu_solver_type


  private :: d_ilu_solver_bld, d_ilu_solver_apply, &
       &  d_ilu_solver_free, &
       &  d_ilu_solver_descr,  d_ilu_solver_sizeof, &
       &  d_ilu_solver_default, d_ilu_solver_dmp, &
       &  d_ilu_solver_apply_vect, d_ilu_solver_get_nzeros, &
       &  d_ilu_solver_get_fmt, d_ilu_solver_check, &
       &  d_ilu_solver_get_id, d_ilu_solver_get_wrksize


  interface 
    subroutine amg_d_ilu_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         & trans,work,wv,info,init,initu)
      import :: psb_desc_type, amg_d_ilu_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)             :: desc_data
      class(amg_d_ilu_solver_type), intent(inout) :: sv
      type(psb_d_vect_type),intent(inout)         :: x
      type(psb_d_vect_type),intent(inout)         :: y
      real(psb_dpk_),intent(in)                    :: alpha,beta
      character(len=1),intent(in)                   :: trans
      real(psb_dpk_),target, intent(inout)         :: work(:)
      type(psb_d_vect_type),intent(inout)         :: wv(:)
      integer(psb_ipk_), intent(out)                :: info
      character, intent(in), optional                :: init
      type(psb_d_vect_type),intent(inout), optional   :: initu
    end subroutine amg_d_ilu_solver_apply_vect
  end interface

  interface 
    subroutine amg_d_ilu_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, amg_d_ilu_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type, psb_ipk_
      implicit none 
      type(psb_desc_type), intent(in)      :: desc_data
      class(amg_d_ilu_solver_type), intent(inout) :: sv
      real(psb_dpk_),intent(inout)         :: x(:)
      real(psb_dpk_),intent(inout)         :: y(:)
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)           :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer(psb_ipk_), intent(out)        :: info
      character, intent(in), optional       :: init
      real(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine amg_d_ilu_solver_apply
  end interface

  interface 
    subroutine amg_d_ilu_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, amg_d_ilu_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none 
      type(psb_dspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(amg_d_ilu_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_dspmat_type), intent(in), target, optional :: b
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_d_ilu_solver_bld
  end interface

  interface 
    subroutine amg_d_ilu_solver_cnv(sv,info,amold,vmold,imold)
      import :: amg_d_ilu_solver_type, psb_dpk_, &
           & psb_d_base_sparse_mat, psb_d_base_vect_type,&
           & psb_ipk_, psb_i_base_vect_type
      implicit none 
      class(amg_d_ilu_solver_type), intent(inout)         :: sv
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional   :: imold
    end subroutine amg_d_ilu_solver_cnv
  end interface
  
  interface 
    subroutine amg_d_ilu_solver_dmp(sv,desc,level,info,prefix,head,solver,global_num)
      import :: psb_desc_type, amg_d_ilu_solver_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type, &
           & psb_ipk_
      implicit none 
      class(amg_d_ilu_solver_type), intent(in) :: sv
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(in)              :: level
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional     :: prefix, head
      logical, optional, intent(in)              :: solver, global_num
    end subroutine amg_d_ilu_solver_dmp
  end interface
  
  interface
    subroutine amg_d_ilu_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, &
           & amg_d_base_solver_type, amg_d_ilu_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_d_ilu_solver_type), intent(inout)               :: sv
      class(amg_d_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)              :: info
    end subroutine amg_d_ilu_solver_clone
  end interface

  interface
    subroutine amg_d_ilu_solver_clone_settings(sv,svout,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, &
           & amg_d_base_solver_type, amg_d_ilu_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_d_ilu_solver_type), intent(inout) :: sv
      class(amg_d_base_solver_type), intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine amg_d_ilu_solver_clone_settings
  end interface

  interface
    subroutine amg_d_ilu_solver_clear_data(sv,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, &
           & amg_d_ilu_solver_type, psb_ipk_
      Implicit None
      
      ! Arguments
      class(amg_d_ilu_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine amg_d_ilu_solver_clear_data
  end interface

contains

  subroutine d_ilu_solver_default(sv)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(inout) :: sv

    sv%fact_type = psb_ilu_n_
    sv%fill_in   = 0
    sv%thresh    = dzero

    return
  end subroutine d_ilu_solver_default

  subroutine d_ilu_solver_check(sv,info)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_) :: err_act
    character(len=20) :: name='d_ilu_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call amg_check_def(sv%fact_type,&
         & 'Factorization',psb_ilu_n_,is_legal_ilu_fact)

    select case(sv%fact_type)
    case(psb_ilu_n_,psb_milu_n_)      
      call amg_check_def(sv%fill_in,&
           & 'Level',izero,is_int_non_negative)
    case(psb_ilu_t_)                 
      call amg_check_def(sv%thresh,&
           & 'Eps',dzero,is_legal_d_fact_thrs)
    end select
    
    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine d_ilu_solver_check

  subroutine d_ilu_solver_cseti(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(inout) :: sv 
    character(len=*), intent(in)                  :: what 
    integer(psb_ipk_), intent(in)                 :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_ilu_solver_cseti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(psb_toupper(trim((what))))
    case('SUB_SOLVE') 
      sv%fact_type = val
    case('SUB_FILLIN')
      sv%fill_in   = val
    case default
      call sv%amg_d_base_solver_type%set(what,val,info,idx=idx)
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_ilu_solver_cseti

  subroutine d_ilu_solver_csetc(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(inout) :: sv
    character(len=*), intent(in)                  :: what 
    character(len=*), intent(in)                  :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act, ival
    character(len=20)  :: name='d_ilu_solver_csetc'

    info = psb_success_
    call psb_erractionsave(err_act)
    ival = amg_stringval(val)
    select case(psb_toupper(trim((what))))
    case('SUB_SOLVE') 
      sv%fact_type = ival
    case default
      call sv%amg_d_base_solver_type%set(what,val,info,idx=idx)
    end select
      
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_ilu_solver_csetc
  
  subroutine d_ilu_solver_csetr(sv,what,val,info,idx)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(inout) :: sv 
    character(len=*), intent(in)                  :: what 
    real(psb_dpk_), intent(in)                     :: val
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), intent(in), optional       :: idx
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_ilu_solver_csetr'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case(psb_toupper(what))
    case('SUB_ILUTHRS') 
      sv%thresh = val
    case default
      call sv%amg_d_base_solver_type%set(what,val,info,idx=idx)
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_ilu_solver_csetr

  subroutine d_ilu_solver_free(sv,info)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(inout) :: sv
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='d_ilu_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_

    if (allocated(sv%d)) then 
      deallocate(sv%d,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if
    call sv%l%free()
    call sv%u%free()
    call sv%dv%free(info)

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_ilu_solver_free

  subroutine d_ilu_solver_descr(sv,info,iout,coarse)

    Implicit None

    ! Arguments
    class(amg_d_ilu_solver_type), intent(in) :: sv
    integer(psb_ipk_), intent(out)             :: info
    integer(psb_ipk_), intent(in), optional    :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer(psb_ipk_) :: err_act
    character(len=20), parameter :: name='amg_d_ilu_solver_descr'
    integer(psb_ipk_) :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = psb_out_unit
    endif

    write(iout_,*) '  Incomplete factorization solver: ',&
         &  amg_fact_names(sv%fact_type)
    select case(sv%fact_type)
    case(psb_ilu_n_,psb_milu_n_)      
      write(iout_,*) '  Fill level:',sv%fill_in
    case(psb_ilu_t_)         
      write(iout_,*) '  Fill level:',sv%fill_in
      write(iout_,*) '  Fill threshold :',sv%thresh
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine d_ilu_solver_descr

  function d_ilu_solver_get_nzeros(sv) result(val)

    implicit none 
    ! Arguments
    class(amg_d_ilu_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_) :: i
    
    val = 0 
    val = val + sv%dv%get_nrows()
    val = val + sv%l%get_nzeros()
    val = val + sv%u%get_nzeros()

    return
  end function d_ilu_solver_get_nzeros

  function d_ilu_solver_sizeof(sv) result(val)

    implicit none 
    ! Arguments
    class(amg_d_ilu_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer(psb_ipk_) :: i

    val = 2*psb_sizeof_ip + psb_sizeof_dp
    val = val + sv%dv%sizeof()
    val = val + sv%l%sizeof()
    val = val + sv%u%sizeof()

    return
  end function d_ilu_solver_sizeof

  function d_ilu_solver_get_fmt() result(val)
    implicit none 
    character(len=32)  :: val

    val = "ILU solver"
  end function d_ilu_solver_get_fmt

  function d_ilu_solver_get_id() result(val)
    implicit none 
    integer(psb_ipk_)  :: val
    
    val = psb_ilu_n_
  end function d_ilu_solver_get_id

  function d_ilu_solver_get_wrksize() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = 2
  end function d_ilu_solver_get_wrksize
  
end module amg_d_ilu_solver

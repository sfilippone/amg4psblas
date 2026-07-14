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
!         software without specific prior written permission.
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
! File: amg_cmlprec_aply.f90
!
! Subroutine: amg_cmlprec_aply
! Version:    real
!
!  Current version of this file contributed by:
!        Ambra Abdullahi Hassan 
!
!
!  This routine computes
!  
!                        Y = beta*Y + alpha*op(ML^(-1))*X,
!  where 
!  - ML is a multilevel preconditioner associated with
!    a certain matrix A and stored in p,
!  - op(ML^(-1)) is ML^(-1) or its transpose, according to the value of trans,
!  - X and Y are vectors,
!  - alpha and beta are scalars.
!
!  The following multilevel strategies can be applied:
!
!  - Additive multilevel Schwarz,
!  - classical V-cycle,
!  - classical W-cycle,
!  - K-cycle both for symmetric and nonsymmetric matrices, where 2 iterations
!    of FCG(1) or GCR, respectively, are applied at each level
!    except the coarsest.
!
!  For each level we have as many submatrices as processes (except for the coarsest
!  level where we might have a replicated index space) and each process takes care
!  of one submatrix.
!
!  A multilevel preconditioner is regarded as an array of 'one-level' data structures,
!  each containing the part of the preconditioner associated to a certain level
!  (for more details see the description of amg_Tonelev_type in amg_prec_type.f90).
!  For each level lev, there is a smoother stored in
!     p%precv(lev)%sm
!  which in turn contains a solver
!     p$precv(lev)%sm%sv
!  Typically the solver acts only locally, and the smoother applies any required
!  parallel communication/action. 
!  Each level has  a matrix A(lev), obtained by 'tranferring' the original
!  matrix A (i.e. the matrix to be preconditioned) to the level lev, through smoothed
!  aggregation.
!
!  The levels are numbered in increasing order starting from the finest one, i.e.
!  level 1 is the finest level and A(1) is the matrix A.
!
!  This routine is formulated in a recursive way, so it is quite compact.
!
!  The V-cycle can be described as follows, where
!  P(lev) denotes the smoothed prolongator from level lev to level
!  lev-1, while R(lev) denotes the corresponding  restriction operator
!  (normally its transpose) from level lev-1 to level lev.
!  M(lev) is the smoother at the current level.
!
!
!   1. Transfer the outer vector Xest to u(1) (inner X at level 1)
!
!   2. Invoke V-cycle(1,M,P,R,A,b,u)
!
!    procedure V-cycle(lev,M,P,R,A,b,u)
!
!      if (lev < nlev) then
!
!         u(lev)   = u(lev) + M(lev)*(b(lev)-A(lev)*u(lev))
!
!         b(lev+1) = R(lev+1)*(b(lev)-A(lev)*u(lev))
!
!         u(lev+1) = V-cycle(lev+1,M,P,R,A,b,u)
!
!         u(lev)   = u(lev) + P(lev+1) * u(lev+1)
!
!         u(lev)   = u(lev) + M(lev)*(b(lev)-A(lev)*u(lev))
!
!      else
!
!         solve  A(lev)*u(lev) = b(lev)
!
!      end if
!
!      return u(lev)
!    end
!
!   3. Transfer u(1) to the external:
!         Yext = beta*Yext + alpha*u(1)
!
!
!  In the implementation, the recursive procedure is inner_ml_aply, which
!  in turn uses amg_inner_add (for additive multilevel),
!  amg_inner_mult (for V-cycle and W-cycle), and
!  amg_inner_k_cycle (for symmetric and non-symmetric K-cycle).
!  
!  For a detailed description of the algorithms, see:
!
!  - B.F. Smith, P.E. Bjorstad, W.D. Gropp,
!    Domain decomposition: parallel multilevel methods for elliptic partial
!    differential equations, Cambridge University Press, 1996.
!
!  - W. L. Briggs, V. E. Henson, S. F.  McCormick,
!    A Multigrid Tutorial, Second Edition
!    SIAM, 2000.
!
!  - K. Stuben,
!    An Introduction to Algebraic Multigrid,
!    in A. Schuller, U. Trottenberg, C. Oosterlee, Multigrid, Academic Press, 2001.
!
!  - Y. Notay, P. S. Vassilevski,
!    Recursive Krylov-based multigrid cycles
!    Numerical Linear Algebra with Applications, 15 (5), 2008, 473--487.
!
!
! Arguments:
!   alpha      -   complex(psb_spk_), input.
!                  The scalar alpha.
!   p          -   type(amg_cprec_type), input.
!                  The multilevel preconditioner data structure containing the
!                  local part of the preconditioner to be applied.
!      Note that nlev = size(p%precv) = number of levels.
!      p%precv(lev)%sm        -  type(psb_cbaseprec_type)
!                                 The pre-'smoother' for the current level
!      p%precv(lev)%sm2       -  type(psb_cbaseprec_type)
!                                 The post-'smoother' for the current level
!                                 may be the same or different from %sm
!      p%precv(lev)%ac        -  type(psb_cspmat_type) 
!                                 The local part of the matrix A(lev).
!      p%precv(lev)%parms     -  type(psb_sml_parms) 
!                                 Parameters controllin the multilevel prec. 
!      p%precv(lev)%desc_ac   -  type(psb_desc_type).
!                                 The communication descriptor associated to the sparse
!                                 matrix A(lev)
!      p%precv(lev)%map       -  type(psb_inter_desc_type)
!                                 Stores the linear operators mapping level (lev-1)
!                                 to (lev) and vice versa. These are the restriction
!                                 and prolongation operators described in the sequel. 
!      p%precv(lev)%base_a    -  type(psb_cspmat_type), pointer.
!                                 Pointer (really a pointer!) to the base matrix of
!                                 the current level, i.e. the local part of A(lev);
!                                 so we have a unified treatment of residuals. We
!                                 need this to avoid passing explicitly the matrix
!                                 A(lev) to the routine which applies the
!                                 preconditioner.
!      p%precv(lev)%base_desc -  type(psb_desc_type), pointer.
!                                 Pointer to the communication descriptor associated
!                                 to the sparse matrix pointed by base_a.  
!                  
!   x          -  complex(psb_spk_), dimension(:), input.
!                 The local part of the vector X.
!   beta       -  complex(psb_spk_), input.
!                 The scalar beta.
!   y          -  complex(psb_spk_), dimension(:), input/output.
!                 The local part of the vector Y.
!   desc_data  -  type(psb_desc_type), input.
!                 The communication descriptor associated to the matrix to be
!                 preconditioned.
!   trans      -  character, optional.
!                 If trans='N','n' then op(M^(-1)) = M^(-1);
!                 if trans='T','t' then op(M^(-1)) = M^(-T) (transpose of M^(-1)).
!   work       -  complex(psb_spk_), dimension (:), optional, target.
!                 Workspace. Its size must be at least 4*desc_data%get_local_cols().
!   info       -  integer, output.
!                 Error code.
!
!   Note that when the LU factorization of the matrix A(lev) is computed instead of
!   the ILU one, by using UMFPACK or SuperLU or MUMPS, the corresponding
!   L and U factors are stored in data structures handled
!   by the third party software. 
!

!
! Old routine for arrays instead of psb_X_vector. To be deleted eventually.
!
!
subroutine amg_cmlprec_aply_a(alpha,p,x,beta,y,desc_data,trans,work,info)

  use psb_base_mod
  use amg_base_prec_type
  use amg_c_inner_mod, amg_protect_name => amg_cmlprec_aply_a

  implicit none

  ! Arguments
  type(psb_desc_type),intent(in)    :: desc_data
  type(amg_cprec_type), intent(inout)  :: p
  complex(psb_spk_),intent(in)         :: alpha,beta
  complex(psb_spk_),intent(inout)      :: x(:)
  complex(psb_spk_),intent(inout)      :: y(:)
  character, intent(in)             :: trans
  complex(psb_spk_),target             :: work(:)
  integer(psb_ipk_), intent(out)              :: info

  ! Local variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me
  integer(psb_ipk_)   :: err_act
  integer(psb_ipk_)   :: debug_level, debug_unit, nlev,nc2l,nr2l,level
  character(len=20)   :: name
  character           :: trans_
  type amg_mlwrk_type
    complex(psb_spk_), allocatable  :: tx(:), ty(:), x2l(:), y2l(:)
  end type amg_mlwrk_type
  type(amg_mlwrk_type), allocatable, target  :: mlwrk(:)

  name='amg_cmlprec_aply'
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Entry  ', size(p%precv)

  trans_ = psb_toupper(trans)

  nlev = size(p%precv)
  allocate(mlwrk(nlev),stat=info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  level = 1

  do level = 1, nlev
    call psb_geasb(mlwrk(level)%x2l,&
         & p%precv(level)%base_desc,info)
    call psb_geasb(mlwrk(level)%y2l,&
         & p%precv(level)%base_desc,info)
    call psb_geasb(mlwrk(level)%tx,&
         & p%precv(level)%base_desc,info)
    call psb_geasb(mlwrk(level)%ty,&
         & p%precv(level)%base_desc,info)
    if (psb_errstatus_fatal()) then 
      nc2l = p%precv(level)%base_desc%get_local_cols()
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
           & a_err='complex(psb_spk_)')
      goto 9999      
    end if
  end do

  mlwrk(level)%x2l(:) = x(:) 
  mlwrk(level)%y2l(:) = czero 

  call inner_ml_aply(level,p,mlwrk,trans_,work,info)    

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Inner prec aply')
    goto 9999
  end if

  call psb_geaxpby(alpha,mlwrk(level)%y2l,beta,y,&
       &   p%precv(level)%base_desc,info)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,&
         & a_err='Error final update')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

    !
  !
  ! inner_ml_aply: apply AMG at a given level.
  ! This routine dispatches the computation according to the type
  ! specified at the current level.
  ! Each of the corrections will inturn call recursively this routine.
  !
  ! Assumptions:
  ! On input:
  !   mlprec_wkr(level)%vx2l   contains the input vector (RHS)
  !   mlprec_wkr(level)%vy2l   contains the initial guess
  !
  ! On output:
  !   mlprec_wkr(level)%vy2l   contains the solution
  !
  ! Constraints: each of the called routines must properly handle
  ! the input/output conditions for level+1 (i.e. apply
  ! prolongation/restriction).
  ! Note: for historical/convenience reasons the prolongator/restrictor
  ! between level and level+1 are stored at level+1. 
  !
  !
  recursive subroutine inner_ml_aply(level,p,mlwrk,trans,work,info)    

    implicit none 

    ! Arguments
    integer(psb_ipk_)                           :: level 
    type(amg_cprec_type), target, intent(inout) :: p
    type(amg_mlwrk_type), intent(inout), target    :: mlwrk(:)
    character, intent(in)                       :: trans
    complex(psb_spk_),target                      :: work(:)
    integer(psb_ipk_), intent(out)              :: info

    type(psb_c_vect_type) :: res
    type(psb_c_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: np, me
    integer(psb_ipk_)   :: i, err_act
    integer(psb_ipk_)   :: debug_level, debug_unit
    integer(psb_ipk_)   :: nlev, ilev, sweeps
    logical             :: pre, post
    character(len=20)   :: name



    name = 'inner_ml_aply'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_ml')
      goto 9999      
    end if
    ctxt = p%precv(level)%base_desc%get_context()
    call psb_info(ctxt, me, np)

    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_ml_aply at level ',level
    end if

    select case(p%precv(level)%parms%ml_cycle) 

    case(amg_no_ml_)
      !
      ! No preconditioning, should not really get here
      ! 
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='amg_no_ml_ in mlprc_aply?')
      goto 9999      

    case(amg_add_ml_)

      call amg_c_inner_add(p, mlwrk, level, trans, work)

    case(amg_mult_ml_, amg_vcycle_ml_, amg_wcycle_ml_)

      call amg_c_inner_mult(p, mlwrk, level, trans, work)
      
! !$    case(amg_kcycle_ml_, amg_kcyclesym_ml_)
! !$
! !$      call amg_c_inner_k_cycle(p, mlwrk, level, trans, work)
      
    case default
      info = psb_err_from_subroutine_ai_
      call psb_errpush(info,name,a_err='invalid ml_cycle',&
           &  i_Err=(/p%precv(level)%parms%ml_cycle,izero,izero,izero,izero/))
      goto 9999      

    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine inner_ml_aply

  recursive subroutine amg_c_inner_add(p, mlwrk, level, trans, work)
    use psb_base_mod
    use amg_prec_mod

    implicit none

    !Input/Oputput variables
    type(amg_cprec_type), intent(inout)  :: p

    type(amg_mlwrk_type), target,  intent(inout) :: mlwrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans
    complex(psb_spk_),target            :: work(:)
    type(psb_c_vect_type) :: res
    type(psb_c_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: np, me
    integer(psb_ipk_)   :: i, err_act
    integer(psb_ipk_)   :: debug_level, debug_unit
    integer(psb_ipk_)   :: nlev, ilev, sweeps
    logical             :: pre, post
    character(len=20)   :: name



    name = 'inner_inner_add'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_add')
      goto 9999      
    end if
    ctxt = p%precv(level)%base_desc%get_context()
    call psb_info(ctxt, me, np)

    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_add at level ',level
    end if

    if ((level<1).or.(level>nlev)) then
      info = psb_err_internal_error_ 
      call psb_errpush(info,name,&
           & a_err='Invalid LEVEL>NLEV')
      goto 9999
    end if
    
    sweeps = p%precv(level)%parms%sweeps_pre
    call p%precv(level)%sm%apply(cone,&
         & mlwrk(level)%x2l,czero,mlwrk(level)%y2l,&
         & p%precv(level)%base_desc, trans,&
         & sweeps,work,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during ADD smoother_apply')
      goto 9999
    end if

    if (level < nlev) then
      ! Apply the restriction
      call p%precv(level+1)%map_rstr(cone,mlwrk(level)%x2l,&
           & czero,mlwrk(level+1)%x2l,&
           & info,work=work)
      mlwrk(level+1)%y2l(:) = czero
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during restriction')
        goto 9999
      end if
      
      call inner_ml_aply(level+1,p,mlwrk,trans,work,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in recursive call')
        goto 9999
      end if

      !
      ! Apply the prolongator and add correction.
      !  
      call p%precv(level+1)%map_prol(cone,&
           & mlwrk(level+1)%y2l,cone,mlwrk(level)%y2l,&
           & info,work=work)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during prolongation')
        goto 9999
      end if


    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_c_inner_add

  recursive subroutine amg_c_inner_mult(p, mlwrk, level, trans, work)
    use psb_base_mod
    use amg_prec_mod

    implicit none

    !Input/Oputput variables
    type(amg_cprec_type), intent(inout)  :: p

    type(amg_mlwrk_type), target,  intent(inout) :: mlwrk(:)
    integer(psb_ipk_), intent(in) :: level
    character, intent(in)             :: trans
    complex(psb_spk_),target            :: work(:)
    type(psb_c_vect_type) :: res
    type(psb_c_vect_type), pointer :: current
    integer(psb_ipk_) :: sweeps_post, sweeps_pre
    ! Local variables
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: np, me
    integer(psb_ipk_)   :: i, err_act
    integer(psb_ipk_)   :: debug_level, debug_unit
    integer(psb_ipk_)   :: nlev, ilev, sweeps
    logical             :: pre, post
    character(len=20)   :: name



    name = 'inner_inner_mult'
    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    nlev = size(p%precv)
    if ((level < 1) .or. (level > nlev)) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='wrong call level to inner_mult')
      goto 9999      
    end if
    ctxt = p%precv(level)%base_desc%get_context()
    call psb_info(ctxt, me, np)

    if(debug_level > 1) then
      write(debug_unit,*) me,' inner_mult at level ',level
    end if

    if ((level < nlev).or.(nlev == 1)) then
      sweeps_post = p%precv(level)%parms%sweeps_post
      sweeps_pre  = p%precv(level)%parms%sweeps_pre
    else
      sweeps_post = p%precv(level-1)%parms%sweeps_post
      sweeps_pre  = p%precv(level-1)%parms%sweeps_pre
    endif

    pre  = ((sweeps_pre>0).and.(trans=='N')).or.((sweeps_post>0).and.(trans/='N'))
    post = ((sweeps_post>0).and.(trans=='N')).or.((sweeps_pre>0).and.(trans/='N'))

    
    if (level < nlev) then 

      !
      ! Apply the first smoother
      !

      if (pre) then
        if (trans == 'N') then 
          sweeps = p%precv(level)%parms%sweeps_pre
          if (info == psb_success_) call p%precv(level)%sm%apply(cone,&
               & mlwrk(level)%x2l,czero,mlwrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info,init='Y')
        else
          sweeps = p%precv(level)%parms%sweeps_post
          if (info == psb_success_) call p%precv(level)%sm2%apply(cone,&
               & mlwrk(level)%x2l,czero,mlwrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info,init='Y')
        end if
        
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during PRE smoother_apply')
          goto 9999
        end if
      endif

      !
      ! Compute the residual and call recursively
      !
      if (pre) then
        call psb_geaxpby(cone,mlwrk(level)%x2l,&
             & czero,mlwrk(level)%ty,&
             & p%precv(level)%base_desc,info)
        
        if (info == psb_success_) call psb_spmm(-cone,p%precv(level)%base_a,&
             & mlwrk(level)%y2l,cone,mlwrk(level)%ty,&
             & p%precv(level)%base_desc,info,work=work,trans=trans)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during residue')
          goto 9999
        end if
        call p%precv(level+1)%map_rstr(cone,mlwrk(level)%ty,&
             & czero,mlwrk(level+1)%x2l,info,work=work)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during restriction')
          goto 9999
        end if
      else
        ! Shortcut: just transfer x2l. 
        call p%precv(level+1)%map_rstr(cone,mlwrk(level)%x2l,&
             & czero,mlwrk(level+1)%x2l,info,work=work)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during restriction')
          goto 9999
        end if
      endif
      ! First guess is zero
      mlwrk(level+1)%y2l(:) = czero
      
      
      call inner_ml_aply(level+1,p,mlwrk,trans,work,info)
      
      if (p%precv(level)%parms%ml_cycle == amg_wcycle_ml_) then
        ! On second call will use output y2l as initial guess
        if (info == psb_success_) call inner_ml_aply(level+1,p,mlwrk,trans,work,info)
      endif
      
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error in recursive call')
        goto 9999
      end if
      
      
      !
      ! Apply the prolongator
      !  
      call p%precv(level+1)%map_prol(cone,mlwrk(level+1)%y2l,&
           & cone,mlwrk(level)%y2l,info,work=work)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during prolongation')
        goto 9999
      end if
      
      !
      ! Compute the residual
      !
      if (post) then
        call psb_geaxpby(cone,mlwrk(level)%x2l,&
             & czero,mlwrk(level)%tx,&
             & p%precv(level)%base_desc,info)
        call psb_spmm(-cone,p%precv(level)%base_a,mlwrk(level)%y2l,&
             & cone,mlwrk(level)%tx,p%precv(level)%base_desc,info,&
             & work=work,trans=trans)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during residue')
          goto 9999
        end if
        !
        ! Apply the second smoother
        !
        if (trans == 'N') then
          sweeps = p%precv(level)%parms%sweeps_post
          if (info == psb_success_) call p%precv(level)%sm2%apply(cone,&
               & mlwrk(level)%tx,cone,mlwrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info,init='Z')
        else 
          sweeps = p%precv(level)%parms%sweeps_pre
          if (info == psb_success_) call p%precv(level)%sm%apply(cone,&
               & mlwrk(level)%tx,cone,mlwrk(level)%y2l,&
               & p%precv(level)%base_desc, trans,&
               & sweeps,work,info,init='Z')
        end if
        
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_,name,&
               & a_err='Error during POST smoother_apply')
          goto 9999
        end if
        
      endif            
      
    else if (level == nlev) then
      
      sweeps = p%precv(level)%parms%sweeps_pre
      if (info == psb_success_) call p%precv(level)%sm%apply(cone,&
           & mlwrk(level)%x2l,czero,mlwrk(level)%y2l,&
           & p%precv(level)%base_desc, trans,&
           & sweeps,work,info)

    else

      info = psb_err_internal_error_ 
      call psb_errpush(info,name,&
           & a_err='Invalid LEVEL vs NLEV')
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine amg_c_inner_mult
    

end subroutine amg_cmlprec_aply_a

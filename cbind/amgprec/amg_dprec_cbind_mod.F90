module amg_dprec_cbind_mod

  use iso_c_binding
  use amg_prec_mod
  use psb_base_cbind_mod
  use psb_dlinsolve_cbind_mod

  type, bind(c) :: amg_c_dprec
    type(c_ptr) :: item = c_null_ptr
  end type amg_c_dprec

contains

#if 1
#define AMGC_DEBUG(MSG) write(*,*) __FILE__,':',__LINE__,':',MSG
#define AMGC_ERROR(MSG) write(*,*) __FILE__,':',__LINE__,':'," ERROR: ",MSG
#else
#define AMGC_DEBUG(MSG)
#define AMGC_ERROR(MSG)
#endif
#define amg_success_ 0
  !#define AMGC_ERR_FILTER(INFO) min(0,INFO)
#define AMGC_ERR_FILTER(INFO) (INFO)
#define AMGC_ERR_HANDLE(INFO) if(INFO/=amg_success_)AMGC_ERROR("ERROR!")
  function  amg_c_dprecinit(cctxt,ph,ptype) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)  :: res
    type(amg_c_dprec)    :: ph
    type(psb_c_object_type), value :: cctxt
    character(c_char)     :: ptype(*)
    integer(psb_ipk_)     :: iret
    type(amg_dprec_type), pointer :: precp
    character(len=80)     :: fptype

    res = -1
    if (c_associated(ph%item)) then
      res = 0
      return
    end if

    allocate(precp,stat=iret)
    if (iret /= 0) return

    ph%item = c_loc(precp)

    call stringc2f(ptype,fptype)

    call precp%init(psb_c2f_ctxt(cctxt),fptype,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecinit

  function  amg_c_dprecseti(ph,what,val) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    integer(psb_c_ipk_), value :: val
    integer(psb_ipk_)     :: iret
    character(len=80)     :: fwhat
    type(amg_dprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)

    call precp%set(fwhat,val,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecseti


  function  amg_c_dprecsetr(ph,what,val) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    real(c_double), value :: val
    integer(psb_ipk_)     :: iret
    character(len=80)     :: fwhat
    type(amg_dprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)

    call precp%set(fwhat,val,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecsetr

  function  amg_c_dprecsetc(ph,what,val) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*), val(*)
    integer(psb_ipk_)     :: iret
    character(len=80)     :: fwhat,fval
    type(amg_dprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)
    call stringc2f(val,fval)

    call precp%set(fwhat,fval,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecsetc

  function  amg_c_dprecbld(ah,cdh,ph) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(len=80)     :: fptype
    integer(psb_ipk_)     :: iret

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call amg_precbld(ap,descp,precp,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)

    return
  end function amg_c_dprecbld

  function  amg_c_dhierarchy_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(len=80)     :: fptype
    integer(psb_ipk_)     :: iret, act


    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call precp%hierarchy_build(ap,descp,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    if (res /=0) then
      act = psb_act_abort_
      call psb_error_handler(act)
    end if

    return
  end function amg_c_dhierarchy_build

  function  amg_c_dsmoothers_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(len=80)     :: fptype
    integer(psb_ipk_)     :: iret, act

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call precp%smoothers_build(ap,descp,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    if (res /=0) then
      act = psb_act_abort_
      call psb_error_handler(act)
    end if

    return
  end function amg_c_dsmoothers_build

  function  amg_c_dkrylov(methd,&
       & ah,ph,bh,xh,cdh,s1,s2,options) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_linsolve_mod
    use psb_prec_cbind_mod
    use psb_dlinsolve_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh,s1,s2
    character(c_char)       :: methd(*)
    type(solveroptions)     :: options


    res= amg_c_dkrylov_opt(methd, ah, ph, bh, xh, options%eps,cdh,  &
         & itmax=options%itmax, iter=options%iter,&
         & itrace=options%itrace, istop=options%istop,&
         & irst=options%irst, err=options%err, s1=s1,s2=s2)

  end function amg_c_dkrylov


  function  amg_c_dkrylov_opt(methd,&
       & ah,ph,bh,xh,eps,cdh,itmax,iter,err,&
       & itrace,irst,istop,s1,s2) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_linsolve_mod
    use psb_objhandle_mod
    use psb_prec_cbind_mod
    use psb_base_string_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh,s1, s2
    integer(psb_c_ipk_), value :: itmax,itrace,irst,istop
    real(c_double), value :: eps
    integer(psb_c_ipk_)        :: iter
    real(c_double)        :: err
    character(c_char)       :: methd(*)

    type(psb_desc_type), pointer   :: descp
    type(psb_dspmat_type), pointer :: ap
    type(amg_dprec_type), pointer  :: precp
    type(psb_d_vect_type), pointer :: xp, bp, s1p, s2p

    integer(psb_ipk_)     :: iret,fitmax,fitrace,first,fistop,fiter
    character(len=20)     :: fmethd
    real(kind(1.d0))      :: feps,ferr

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(bh%item)) then
      call c_f_pointer(bh%item,bp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if
    if (c_associated(s1%item)) then
      call c_f_pointer(s1%item,s1p)
    else
      nullify(s1p)
    end if
    if (c_associated(s2%item)) then
      call c_f_pointer(s2%item,s2p)
    else
      nullify(s2p)
    end if
    write(0,*) 'krylov: s1 s2 ',associated(s1p),associated(s2p)
    call stringc2f(methd,fmethd)
    feps    = eps
    fitmax  = itmax
    fitrace = itrace
    first   = irst
    fistop  = istop
!!$    flush(0)
!!$    
!!$    if (associated(s1p)) write(0,*) 'amg_c_dkrylov_opt S1 ',s1p%get_nrows(), &
!!$         & s1p%nrm2(s1p%get_nrows())
!!$    if (associated(s2p)) write(0,*) 'amg_c_dkrylov_opt S2 ',s2p%get_nrows(), &
!!$         & s2p%nrm2(s2p%get_nrows())
!!$      flush(0)
    if (associated(s1p).and.associated(s2p)) then
      call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
           & descp, iret,&
           & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
           & irst=first, err=ferr,s1=s1p,s2=s2p)
    else if  (associated(s1p)) then
      call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
           & descp, iret,&
           & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
           & irst=first, err=ferr,s1=s1p)
    else  if (associated(s2p)) then
      call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
           & descp, iret,&
           & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
           & irst=first, err=ferr,s2=s2p)
    else
      call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
           & descp, iret,&
           & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
           & irst=first, err=ferr)
    end if
    
!!$    call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
!!$         & descp, iret,&
!!$         & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
!!$         & irst=first, err=ferr)

    iter = fiter
    err  = ferr
    res = min(iret,0)

  end function amg_c_dkrylov_opt

  function amg_c_dprecapply(ph,bc,xc,cdh) bind(c,name="amg_c_dprecapply") result(res)
    use psb_base_mod
    use psb_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph ! C handle to preconditioner
    type(psb_c_object_type) :: bc ! C handle to rhs
    type(psb_c_object_type) :: xc ! C handle to solution
    type(psb_c_object_type) :: cdh ! C handle to descriptor

    ! Fortran containers for preconditioner, lhs, rhs and descriptor
    type(amg_dprec_type), pointer  :: precp
    type(psb_d_vect_type), pointer :: xp, bp
    type(psb_desc_type), pointer   :: descp
    integer(psb_ipk_)     :: info

    res = -1
    ! Check descriptor
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    ! Check rhs and solution
    if (c_associated(bc%item)) then
      call c_f_pointer(bc%item,bp)
    else
      return
    end if
    if (c_associated(xc%item)) then
      call c_f_pointer(xc%item,xp)
    else
      return
    end if
    ! Check preconditioner
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if
    ! Apply preconditioner
    call precp%apply(bp,xp,descp,info)

    ! Error handling and return
    res = AMGC_ERR_FILTER(info)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecapply

  function amg_c_dprecapply_opt(ph,bc,xc,cdh,ctrans) bind(c,name="amg_c_dprecapply_opt") result(res)
    use psb_base_mod
    use psb_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph ! C handle to preconditioner
    type(psb_c_object_type) :: bc ! C handle to rhs
    type(psb_c_object_type) :: xc ! C handle to solution
    type(psb_c_object_type) :: cdh ! C handle to descriptor
    character(c_char)       :: ctrans(*) ! Tranpose flag as character

    ! Fortran containers for preconditioner, lhs, rhs and descriptor
    type(amg_dprec_type), pointer  :: precp
    type(psb_d_vect_type), pointer :: xp, bp
    type(psb_desc_type), pointer   :: descp
    character(len=10)      :: ftrans
    integer(psb_ipk_)     :: info

    res = -1
    ! Check descriptor
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    ! Check rhs and solution
    if (c_associated(bc%item)) then
      call c_f_pointer(bc%item,bp)
    else
      return
    end if
    if (c_associated(xc%item)) then
      call c_f_pointer(xc%item,xp)
    else
      return
    end if
    ! Check preconditioner
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    ! Convert transpose flag
    call stringc2f(ctrans,ftrans)

    ! Apply preconditioner
    call precp%apply(bp,xp,descp,info,trans=ftrans)

    ! Error handling and return
    res = AMGC_ERR_FILTER(info)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecapply_opt

  function  amg_c_dprecfree(ph) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    integer(psb_ipk_)     :: iret
    type(amg_dprec_type), pointer :: precp
    character(len=80)     :: fptype

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if


    call precp%free(iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_dprecfree

  function amg_c_ddescr(ph) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    integer(psb_c_ipk_)     :: iret
    type(amg_dprec_type), pointer :: precp

    res = -1
    iret = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if


    call precp%descr(iret)
    call flush(psb_out_unit)

    iret = 0
    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_ddescr

end module amg_dprec_cbind_mod

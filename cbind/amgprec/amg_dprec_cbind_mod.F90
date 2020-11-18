module amg_dprec_cbind_mod

  use iso_c_binding
  use amg_prec_mod
  use psb_base_cbind_mod

  type, bind(c) :: amg_c_dprec
    type(c_ptr) :: item = c_null_ptr
  end type amg_c_dprec

contains

#if 1
#define MLDC_DEBUG(MSG) write(*,*) __FILE__,':',__LINE__,':',MSG
#define MLDC_ERROR(MSG) write(*,*) __FILE__,':',__LINE__,':'," ERROR: ",MSG
#else
#define MLDC_DEBUG(MSG)
#define MLDC_ERROR(MSG)
#endif
#define amg_success_ 0
!#define MLDC_ERR_FILTER(INFO) min(0,INFO)
#define MLDC_ERR_FILTER(INFO) (INFO)
#define MLDC_ERR_HANDLE(INFO) if(INFO/=amg_success_)MLDC_ERROR("ERROR!")

  function  amg_c_dprecinit(cctxt,ph,ptype) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_)  :: res
    type(amg_c_dprec)    :: ph
    type(psb_c_object_type), value :: cctxt
    character(c_char)     :: ptype(*)
    integer               :: info
    type(amg_dprec_type), pointer :: precp
    character(len=80)     :: fptype

    res = -1
    if (c_associated(ph%item)) then
      res = 0
      return
    end if

    allocate(precp,stat=info)
    if (info /= 0) return

    ph%item = c_loc(precp)

    call stringc2f(ptype,fptype)

    call precp%init(psb_c2f_ctxt(cctxt),fptype,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function amg_c_dprecinit

  function  amg_c_dprecseti(ph,what,val) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    integer(psb_c_ipk_), value :: val
    integer               :: info
    character(len=80)     :: fwhat
    type(amg_dprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)

    call amg_precset(precp,fwhat,val,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function amg_c_dprecseti


  function  amg_c_dprecsetr(ph,what,val) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    real(c_double), value :: val
    integer               :: info
    character(len=80)     :: fwhat
    type(amg_dprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)

    call amg_precset(precp,fwhat,val,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function amg_c_dprecsetr

  function  amg_c_dprecsetc(ph,what,val) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*), val(*)
    integer               :: info
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

    call amg_precset(precp,fwhat,fval,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function amg_c_dprecsetc

  function  amg_c_dprecbld(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer               :: info
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(len=80)     :: fptype

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

    call amg_precbld(ap,descp,precp,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dprecbld

  function  amg_c_dhierarchy_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer               :: info
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(len=80)     :: fptype

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

    call precp%hierarchy_build(ap,descp,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dhierarchy_build

  function  amg_c_dsmoothers_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer               :: info
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(len=80)     :: fptype

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

    call precp%smoothers_build(ap,descp,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dsmoothers_build

  function  amg_c_dkrylov(methd,&
       & ah,ph,bh,xh,cdh,options) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_krylov_mod
    use psb_prec_cbind_mod
    use psb_dkrylov_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh
    character(c_char)       :: methd(*)
    type(solveroptions)     :: options

    res= amg_c_dkrylov_opt(methd, ah, ph, bh, xh, options%eps,cdh,  &
         & itmax=options%itmax, iter=options%iter,&
         & itrace=options%itrace, istop=options%istop,&
         & irst=options%irst, err=options%err)

  end function amg_c_dkrylov


  function  amg_c_dkrylov_opt(methd,&
       & ah,ph,bh,xh,eps,cdh,itmax,iter,err,itrace,irst,istop) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_krylov_mod
    use psb_objhandle_mod
    use psb_prec_cbind_mod
    use psb_base_string_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh
    integer(psb_c_ipk_), value :: itmax,itrace,irst,istop
    real(c_double), value :: eps
    integer(psb_c_ipk_)        :: iter
    real(c_double)        :: err
    character(c_char)       :: methd(*)
    type(psb_desc_type), pointer   :: descp
    type(psb_dspmat_type), pointer :: ap
    type(amg_dprec_type), pointer  :: precp
    type(psb_d_vect_type), pointer :: xp, bp

    integer               :: info,fitmax,fitrace,first,fistop,fiter
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


    call stringc2f(methd,fmethd)
    feps    = eps
    fitmax  = itmax
    fitrace = itrace
    first   = irst
    fistop  = istop

    call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
         & descp, info,&
         & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
         & irst=first, err=ferr)
    iter = fiter
    err  = ferr
    res = min(info,0)

  end function amg_c_dkrylov_opt

  function  amg_c_dprecfree(ph) bind(c) result(res)
    use psb_base_mod
    use amg_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    integer               :: info
    type(amg_dprec_type), pointer :: precp
    character(len=80)     :: fptype

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if


    call precp%free(info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function amg_c_dprecfree

  function amg_c_ddescr(ph) bind(c) result(res)
   use psb_base_mod
   use amg_prec_mod
   implicit none
 
   integer(psb_c_ipk_) :: res
   type(psb_c_object_type) :: ph
   integer               :: info
   type(amg_dprec_type), pointer :: precp
 
   res = -1
   info = -1
   if (c_associated(ph%item)) then
     call c_f_pointer(ph%item,precp)
   else
     return
   end if
 
 
   call precp%descr()
   call flush(psb_out_unit)
 
   info = 0
   res = MLDC_ERR_FILTER(info)
   MLDC_ERR_HANDLE(res)
   return
 end function amg_c_ddescr

end module amg_dprec_cbind_mod

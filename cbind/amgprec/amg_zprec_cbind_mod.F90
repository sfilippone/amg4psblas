module mld_zprec_cbind_mod

  use iso_c_binding
  use mld_prec_mod
  use psb_base_cbind_mod

  type, bind(c) :: mld_c_zprec
    type(c_ptr) :: item = c_null_ptr
  end type mld_c_zprec

contains

#if 1
#define MLDC_DEBUG(MSG) write(*,*) __FILE__,':',__LINE__,':',MSG
#define MLDC_ERROR(MSG) write(*,*) __FILE__,':',__LINE__,':'," ERROR: ",MSG
#else
#define MLDC_DEBUG(MSG)
#define MLDC_ERROR(MSG)
#endif
#define mld_success_ 0
!#define MLDC_ERR_FILTER(INFO) min(0,INFO)
#define MLDC_ERR_FILTER(INFO) (INFO)
#define MLDC_ERR_HANDLE(INFO) if(INFO/=mld_success_)MLDC_ERROR("ERROR!")

  function  mld_c_zprecinit(ictxt,ph,ptype) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_)  :: res
    type(mld_c_zprec)    :: ph
    integer(psb_c_ipk_), value :: ictxt
    character(c_char)     :: ptype(*)
    integer               :: info
    type(mld_zprec_type), pointer :: precp
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

    call precp%init(ictxt,fptype,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function mld_c_zprecinit

  function  mld_c_zprecseti(ph,what,val) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    integer(psb_c_ipk_), value :: val
    integer               :: info
    character(len=80)     :: fwhat
    type(mld_zprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)

    call mld_precset(precp,fwhat,val,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function mld_c_zprecseti


  function  mld_c_zprecsetr(ph,what,val) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    real(c_double), value :: val
    integer               :: info
    character(len=80)     :: fwhat
    type(mld_zprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)

    call mld_precset(precp,fwhat,val,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function mld_c_zprecsetr

  function  mld_c_zprecsetc(ph,what,val) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*), val(*)
    integer               :: info
    character(len=80)     :: fwhat,fval
    type(mld_zprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call stringc2f(what,fwhat)
    call stringc2f(val,fval)

    call mld_precset(precp,fwhat,fval,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
    return
  end function mld_c_zprecsetc

  function  mld_c_zprecbld(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer               :: info
    type(mld_zprec_type), pointer  :: precp
    type(psb_zspmat_type), pointer :: ap
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

    call mld_precbld(ap,descp,precp,info)

    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)

    return
  end function mld_c_zprecbld

  function  mld_c_zhierarchy_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer               :: info
    type(mld_zprec_type), pointer  :: precp
    type(psb_zspmat_type), pointer :: ap
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
  end function mld_c_zhierarchy_build

  function  mld_c_zsmoothers_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer               :: info
    type(mld_zprec_type), pointer  :: precp
    type(psb_zspmat_type), pointer :: ap
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
  end function mld_c_zsmoothers_build

  function  mld_c_zkrylov(methd,&
       & ah,ph,bh,xh,cdh,options) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_krylov_mod
    use psb_prec_cbind_mod
    use psb_zkrylov_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh
    character(c_char)       :: methd(*)
    type(solveroptions)     :: options

    res= mld_c_zkrylov_opt(methd, ah, ph, bh, xh, options%eps,cdh,  &
         & itmax=options%itmax, iter=options%iter,&
         & itrace=options%itrace, istop=options%istop,&
         & irst=options%irst, err=options%err)

  end function mld_c_zkrylov


  function  mld_c_zkrylov_opt(methd,&
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
    type(psb_zspmat_type), pointer :: ap
    type(mld_zprec_type), pointer  :: precp
    type(psb_z_vect_type), pointer :: xp, bp

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

  end function mld_c_zkrylov_opt

  function  mld_c_zprecfree(ph) bind(c) result(res)
    use psb_base_mod
    use mld_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    integer               :: info
    type(mld_zprec_type), pointer :: precp
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
  end function mld_c_zprecfree

  function mld_c_zdescr(ph) bind(c) result(res)
   use psb_base_mod
   use mld_prec_mod
   implicit none
 
   integer(psb_c_ipk_) :: res
   type(psb_c_object_type) :: ph
   integer               :: info
   type(mld_zprec_type), pointer :: precp
 
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
 end function mld_c_zdescr

end module mld_zprec_cbind_mod

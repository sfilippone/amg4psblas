module amg_zprec_cbind_mod

  use iso_c_binding
  use amg_prec_mod
  use psb_base_cbind_mod

  type, bind(c) :: amg_c_zprec
    type(c_ptr) :: item = c_null_ptr
  end type amg_c_zprec

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

  function  amg_c_zprecinit(cctxt,ph,ptype) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)  :: res
    type(amg_c_zprec)    :: ph
    type(psb_c_object_type), value :: cctxt
    character(c_char)     :: ptype(*)
    integer(psb_ipk_)     :: iret
    type(amg_zprec_type), pointer :: precp
    character(len=80)     :: fptype

    res = -1
    if (c_associated(ph%item)) then
      res = 0
      return
    end if

    allocate(precp,stat=iret)
    if (iret /= 0) return

    ph%item = c_loc(precp)

    call psb_stringc2f(ptype,fptype)

    call precp%init(psb_c2f_ctxt(cctxt),fptype,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_zprecinit

  function  amg_c_zprecseti(ph,what,val) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    integer(psb_c_ipk_), value :: val
    integer(psb_ipk_)     :: iret
    character(len=80)     :: fwhat
    type(amg_zprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call psb_stringc2f(what,fwhat)

    call precp%set(fwhat,val,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_zprecseti


  function  amg_c_zprecsetr(ph,what,val) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*)
    real(c_double), value :: val
    integer(psb_ipk_)     :: iret
    character(len=80)     :: fwhat
    type(amg_zprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call psb_stringc2f(what,fwhat)

    call precp%set(fwhat,val,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_zprecsetr

  function  amg_c_zprecsetc(ph,what,val) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: what(*), val(*)
    integer(psb_ipk_)     :: iret
    character(len=80)     :: fwhat,fval
    type(amg_zprec_type), pointer  :: precp

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call psb_stringc2f(what,fwhat)
    call psb_stringc2f(val,fval)

    call precp%set(fwhat,fval,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_zprecsetc

  function  amg_c_zprecbld(ah,cdh,ph) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer(psb_ipk_)     :: iret
    type(amg_zprec_type), pointer  :: precp
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

    call amg_precbld(ap,descp,precp,iret)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)

    return
  end function amg_c_zprecbld

  function  amg_c_zhierarchy_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_zprec_type), pointer  :: precp
    type(psb_zspmat_type), pointer :: ap
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
  end function amg_c_zhierarchy_build

  function  amg_c_zsmoothers_build(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_zprec_type), pointer  :: precp
    type(psb_zspmat_type), pointer :: ap
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
  end function amg_c_zsmoothers_build

  function  amg_c_zsmoothers_build_opt(ah,cdh,ph,afmt,cdfmt) bind(c) result(res)
#if defined (PSB_HAVE_CUDA)
    use psb_cuda_mod
#endif
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_zprec_type), pointer  :: precp
    type(psb_zspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(c_char)     :: afmt(*), cdfmt(*)
    character(len=80)     :: fptype
    integer(psb_ipk_)     :: iret
    ! Local variables for formats
    character(len=10)     :: fafmt, fcdfmt
#if defined (PSB_HAVE_CUDA)
    type(psb_z_vect_cuda), target         :: dvgpu
    type(psb_i_vect_cuda), target         :: ivgpu
    ! GPU matrix molds
    type(psb_z_cuda_hlg_sparse_mat), target   :: ahlg
    type(psb_z_cuda_csrg_sparse_mat), target   :: acsrg
    type(psb_z_cuda_elg_sparse_mat), target    :: aelg
#endif
    type(psb_z_base_vect_type), target    :: dvhost
    type(psb_i_base_vect_type), target    :: ivhost
    ! CPU matrix molds
    type(psb_z_ell_sparse_mat), target    :: aell
    type(psb_z_csr_sparse_mat), target   :: acsr
    type(psb_z_coo_sparse_mat), target   :: acoo
    type(psb_z_hll_sparse_mat), target    :: ahll
    type(psb_z_dns_sparse_mat), target   :: adns

    ! molding variables
    class(psb_z_base_vect_type), pointer  :: vmold
    class(psb_z_base_sparse_mat), pointer :: amold
    class(psb_i_base_vect_type), pointer  :: imold


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
    ! Convert formats
    call psb_stringc2f(afmt,fafmt)
    call psb_stringc2f(cdfmt,fcdfmt)
    ! Select matrix mold
    select case (psb_toupper(fafmt))
#if defined (PSB_HAVE_CUDA)
    case('CSRG')
      amold => acsrg
    case('ELG')
      amold => aelg
    case('HLG')
      amold => ahlg
#endif
    case('CSR')
      amold => acsr
    case('ELL')
      amold => aell
    case('COO')
      amold => acoo
    case('HLL')
      amold => ahll
    case('DNS')
      amold => adns
    case default
      write(psb_err_unit,'(A)') 'amg_c_zsmoothers_build_format: Unknown format ', fafmt, ' defaulting to CSR'
      amold => acsr
    end select
    ! Select vector mold
    select case (psb_toupper(fcdfmt))
#if defined (PSB_HAVE_CUDA)
    case('GPU','DEVICE')
      vmold => dvgpu
      imold => ivgpu
#endif
    case('HOST','CPU')
      vmold => dvhost
      imold => ivhost
    case default
      write(psb_err_unit,'(A)') 'amg_c_zsmoothers_build_format: Unknown format ', fcdfmt, ' defaulting to HOST/CPU'
      vmold => dvhost
      imold => ivhost
    end select


    call precp%smoothers_build(ap,descp,iret,amold=amold,vmold=vmold,imold=imold)

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)

    return
  end function amg_c_zsmoothers_build_opt

  function  amg_c_zkrylov(methd,&
       & ah,ph,bh,xh,cdh,s1,s2,options) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_linsolve_mod
    use psb_prec_cbind_mod
    use psb_zlinsolve_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh,s1,s2
    character(c_char)       :: methd(*)
    type(solveroptions)     :: options

    res= amg_c_zkrylov_opt(methd, ah, ph, bh, xh, options%eps,cdh,  &
         & itmax=options%itmax, iter=options%iter,&
         & itrace=options%itrace, istop=options%istop,&
         & irst=options%irst, err=options%err, s1=s1,s2=s2)

  end function amg_c_zkrylov


  function  amg_c_zkrylov_opt(methd,&
       & ah,ph,bh,xh,eps,cdh,itmax,iter,err,&
       & itrace,irst,istop,s1,s2) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_linsolve_mod
    use psb_objhandle_mod
    use psb_prec_cbind_mod
    implicit none
    integer(psb_c_ipk_)          :: res
    type(psb_c_object_type) :: ah,cdh,ph,bh,xh,s1, s2
    integer(psb_c_ipk_), value :: itmax,itrace,irst,istop
    real(c_double), value :: eps
    integer(psb_c_ipk_)        :: iter
    real(c_double)        :: err
    character(c_char)       :: methd(*)
    type(psb_desc_type), pointer   :: descp
    type(psb_zspmat_type), pointer :: ap
    type(amg_zprec_type), pointer  :: precp
    type(psb_z_vect_type), pointer :: xp, bp, s1p, s2p

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


    call psb_stringc2f(methd,fmethd)
    feps    = eps
    fitmax  = itmax
    fitrace = itrace
    first   = irst
    fistop  = istop

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

    iter = fiter
    err  = ferr
    res = min(iret,0)

  end function amg_c_zkrylov_opt

  function amg_c_zprecapply(ph,bc,xc,cdh) bind(c,name="amg_c_zprecapply") result(res)
    use psb_base_mod
    use psb_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph ! C handle to preconditioner
    type(psb_c_object_type) :: bc ! C handle to rhs
    type(psb_c_object_type) :: xc ! C handle to solution
    type(psb_c_object_type) :: cdh ! C handle to descriptor

    ! Fortran containers for preconditioner, lhs, rhs and descriptor
    type(amg_zprec_type), pointer  :: precp
    type(psb_z_vect_type), pointer :: xp, bp
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
  end function amg_c_zprecapply

  function amg_c_zprecapply_opt(ph,bc,xc,cdh,ctrans) bind(c,name="amg_c_zprecapply_opt") result(res)
    use psb_base_mod
    use psb_prec_mod
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph ! C handle to preconditioner
    type(psb_c_object_type) :: bc ! C handle to rhs
    type(psb_c_object_type) :: xc ! C handle to solution
    type(psb_c_object_type) :: cdh ! C handle to descriptor
    character(c_char)       :: ctrans(:) ! Tranpose flag as character

    ! Fortran containers for preconditioner, lhs, rhs and descriptor
    type(amg_zprec_type), pointer  :: precp
    type(psb_z_vect_type), pointer :: xp, bp
    type(psb_desc_type), pointer   :: descp
    character(len=10)     :: ftrans
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
    call psb_stringc2f(ctrans,ftrans)
    
    ! Apply preconditioner
    call precp%apply(bp,xp,descp,info,trans=ftrans)

    ! Error handling and return
    res = AMGC_ERR_FILTER(info)
    AMGC_ERR_HANDLE(res)
    return
  end function amg_c_zprecapply_opt

  function  amg_c_zprecfree(ph) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    integer(psb_ipk_)     :: iret
    type(amg_zprec_type), pointer :: precp
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
  end function amg_c_zprecfree

  function amg_c_zdescr(ph) bind(c) result(res)
   implicit none

   integer(psb_c_ipk_) :: res
   type(psb_c_object_type) :: ph
   integer(psb_c_ipk_)     :: iret
   type(amg_zprec_type), pointer :: precp

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
 end function amg_c_zdescr

   function amg_c_zallocate_wrk(ph,chfmt) bind(c, name="amg_c_zallocate_wrk") result(res)
#if defined (PSB_HAVE_CUDA)
    use psb_cuda_mod
#endif
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: chfmt(*)
    integer(psb_ipk_)     :: iret
    type(amg_zprec_type), pointer :: precp
    character(len=6)     :: fchfmt
! Local variable
    integer(psb_ipk_) :: info
! Local mold variables
#if defined (PSB_HAVE_CUDA)
    type(psb_z_vect_cuda), target         :: zvgpu
#endif
    type(psb_z_base_vect_type), target         :: zvhost
    class(psb_z_base_vect_type), pointer  :: vmold

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call psb_stringc2f(chfmt,fchfmt)
    select case (psb_toupper(fchfmt))
    case('HOST','CPU')
      vmold => zvhost
#if defined (PSB_HAVE_CUDA)
    case('GPU','DEVICE')
      vmold => zvgpu
#endif
    case default
      write(psb_err_unit,'(A)') 'amg_c_zallocate_wrk: Unknown format ', fchfmt, ' defaulting to HOST/CPU'
      vmold => zvhost
    end select

    call precp%allocate_wrk(info,vmold=vmold)

    iret = info

    res = AMGC_ERR_FILTER(iret)
    AMGC_ERR_HANDLE(res)
    return

  end function amg_c_zallocate_wrk

end module amg_zprec_cbind_mod

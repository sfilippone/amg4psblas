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

    call psb_stringc2f(ptype,fptype)

    call precp%init(psb_c2f_ctxt(cctxt),fptype,iret)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
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

    call psb_stringc2f(what,fwhat)

    call precp%set(fwhat,val,iret)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
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

    call psb_stringc2f(what,fwhat)

    call precp%set(fwhat,val,iret)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
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

    call psb_stringc2f(what,fwhat)
    call psb_stringc2f(val,fval)

    call precp%set(fwhat,fval,iret)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
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

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dprecbld

  function  amg_c_dhierarchy_build(ah,cdh,ph) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    integer(psb_ipk_)     :: iret
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

    call precp%hierarchy_build(ap,descp,iret)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dhierarchy_build

  function  amg_c_dsmoothers_build(ah,cdh,ph) bind(c) result(res)
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

    call precp%smoothers_build(ap,descp,iret)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dsmoothers_build

  function  amg_c_dsmoothers_build_format(ah,cdh,ph,afmt,cdfmt) bind(c) result(res)
#if defined (PSB_HAVE_CUDA)
    use psb_cuda_mod
#endif
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type)  :: ph,ah,cdh
    type(amg_dprec_type), pointer  :: precp
    type(psb_dspmat_type), pointer :: ap
    type(psb_desc_type), pointer   :: descp
    character(c_char)     :: afmt(*), cdfmt(*)
    character(len=80)     :: fptype
    integer(psb_ipk_)     :: iret
    ! Local variables for formats
    character(len=10)     :: fafmt, fcdfmt
#if defined (PSB_HAVE_CUDA)
    type(psb_d_vect_cuda), target         :: dvgpu
    type(psb_i_vect_cuda), target         :: ivgpu
    ! GPU matrix molds
    type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
    type(psb_d_cuda_hdiag_sparse_mat), target :: ahdiag
    type(psb_d_cuda_csrg_sparse_mat), target   :: acsrg
    type(psb_d_cuda_elg_sparse_mat), target    :: aelg
#endif
    type(psb_d_base_vect_type), target    :: dvhost
    type(psb_i_base_vect_type), target    :: ivhost
    ! CPU matrix molds
    type(psb_d_ell_sparse_mat), target    :: aell
    type(psb_d_csr_sparse_mat), target   :: acsr
    type(psb_d_coo_sparse_mat), target   :: acoo
    type(psb_d_hll_sparse_mat), target    :: ahll
    type(psb_d_hdia_sparse_mat), target  :: ahdia
    type(psb_d_dns_sparse_mat), target   :: adns

    ! molding variables
    class(psb_d_base_vect_type), pointer  :: vmold
    class(psb_d_base_sparse_mat), pointer :: amold
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
    case('HDIAG')
      amold => ahdiag
#endif
    case('CSR')
      amold => acsr
    case('ELL')
      amold => aell
    case('COO')
      amold => acoo
    case('HLL')
      amold => ahll
    case('HDIA')
      amold => ahdia
    case('DNS')
      amold => adns
    case default
      write(psb_err_unit,'(A)') 'amg_c_dsmoothers_build_format: Unknown format ', fafmt, ' defaulting to CSR'
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
      write(psb_err_unit,'(A)') 'amg_c_dsmoothers_build_format: Unknown format ', fcdfmt, ' defaulting to HOST/CPU'
      vmold => dvhost
      imold => ivhost
    end select


    call precp%smoothers_build(ap,descp,iret,amold=amold,vmold=vmold,imold=imold)

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)

    return
  end function amg_c_dsmoothers_build_format

  function  amg_c_dkrylov(methd,&
       & ah,ph,bh,xh,cdh,options) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_linsolve_mod
    use psb_prec_cbind_mod
    use psb_dlinsolve_cbind_mod
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
    use psb_linsolve_mod
    use psb_objhandle_mod
    use psb_prec_cbind_mod
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


    call psb_stringc2f(methd,fmethd)
    feps    = eps
    fitmax  = itmax
    fitrace = itrace
    first   = irst
    fistop  = istop

    call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
         & descp, iret,&
         & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
         & irst=first, err=ferr)
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
    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
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
    call psb_stringc2f(ctrans,ftrans)
    
    ! Apply preconditioner
    call precp%apply(bp,xp,descp,info,trans=ftrans)

    ! Error handling and return
    res = MLDC_ERR_FILTER(info)
    MLDC_ERR_HANDLE(res)
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

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
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
    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
    return
  end function amg_c_ddescr

  function amg_c_dallocate_wrk(ph,chfmt) bind(c, name="amg_c_dallocate_wrk") result(res)
#if defined (PSB_HAVE_CUDA)
    use psb_cuda_mod
#endif
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: ph
    character(c_char)     :: chfmt(*)
    integer(psb_ipk_)     :: iret
    type(amg_dprec_type), pointer :: precp
    character(len=6)     :: fchfmt
! Local variable
    integer(psb_ipk_) :: info
! Local mold variables
#if defined (PSB_HAVE_CUDA)
    type(psb_d_vect_cuda), target         :: dvgpu
#endif
    type(psb_d_base_vect_type), target         :: dvhost
    class(psb_d_base_vect_type), pointer  :: vmold

    res = -1
    if (c_associated(ph%item)) then
      call c_f_pointer(ph%item,precp)
    else
      return
    end if

    call psb_stringc2f(chfmt,fchfmt)
    select case (psb_toupper(fchfmt))
    case('HOST','CPU')
      vmold => dvhost
#if defined (PSB_HAVE_CUDA)
    case('GPU','DEVICE')
      vmold => dvgpu
#endif
    case default
      write(psb_err_unit,'(A)') 'amg_c_dallocate_wrk: Unknown format ', fchfmt, ' defaulting to HOST/CPU'
      vmold => dvhost
    end select

    call precp%allocate_wrk(info,vmold=vmold)

    iret = info

    res = MLDC_ERR_FILTER(iret)
    MLDC_ERR_HANDLE(res)
    return

  end function amg_c_dallocate_wrk

end module amg_dprec_cbind_mod

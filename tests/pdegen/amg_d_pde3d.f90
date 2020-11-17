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
! File: amg_d_pde3d.f90
!
! Program: amg_d_pde3d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. 
! 
!
! The PDE is a general second order equation in 3d
!
!   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
!      dxdx     dydy       dzdz        dx       dy         dz   
!
! with Dirichlet boundary conditions
!   u = g 
!
!  on the unit cube  0<=x,y,z<=1.
!
!
! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
!
! There are three choices available for data distribution:
! 1. A simple BLOCK distribution
! 2. A ditribution based on arbitrary assignment of indices to processes,
!    typically from a graph partitioner
! 3. A 3D distribution in which the unit cube is partitioned
!    into subcubes, each one assigned to a process.
!
module amg_d_pde3d_mod
  use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_lpk_, psb_desc_type,&
       &  psb_dspmat_type, psb_d_vect_type, dzero,&
       &  psb_d_base_sparse_mat, psb_d_base_vect_type, &
       &  psb_i_base_vect_type, psb_l_base_vect_type

  interface 
    function d_func_3d(x,y,z) result(val)
      import :: psb_dpk_
      real(psb_dpk_), intent(in) :: x,y,z
      real(psb_dpk_) :: val
    end function d_func_3d
  end interface 

  interface amg_gen_pde3d
    module procedure  amg_d_gen_pde3d
  end interface amg_gen_pde3d

  
contains

  function d_null_func_3d(x,y,z) result(val)

    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    
    val = dzero

  end function d_null_func_3d
  !
  ! functions parametrizing the differential equation 
  !  
  !
  ! Note: b1, b2 and b3 are the coefficients of the first
  ! derivative of the unknown function. The default
  ! we apply here is to have them zero, so that the resulting
  ! matrix is symmetric/hermitian and suitable for
  ! testing with CG and FCG.
  ! When testing methods for non-hermitian matrices you can
  ! change the B1/B2/B3 functions to e.g. done/sqrt((3*done))
  !
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=dzero
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=dzero
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      

    b3=dzero
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z      
    c=dzero
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  a1   
    real(psb_dpk_), intent(in) :: x,y,z
    a1=done/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=done/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=done/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then 
      g = exp(y**2-z**2)
    end if
  end function g

  
  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs. 
  !
  subroutine amg_d_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,info,&
       & f,amold,vmold,imold,partition,nrl,iv)
    use psb_base_mod
    use psb_util_mod
    !
    !   Discretizes the partial differential equation
    ! 
    !   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
    ! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
    !      dxdx     dydy       dzdz        dx       dy         dz   
    !
    ! with Dirichlet boundary conditions
    !   u = g 
    !
    !  on the unit cube  0<=x,y,z<=1.
    !
    !
    ! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
    !
    implicit none
    integer(psb_ipk_)     :: idim
    type(psb_dspmat_type) :: a
    type(psb_d_vect_type) :: xv,bv
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type)   :: ctxt
    integer(psb_ipk_)     :: info
    character(len=*)      :: afmt
    procedure(d_func_3d), optional :: f
    class(psb_d_base_sparse_mat), optional :: amold
    class(psb_d_base_vect_type), optional :: vmold 
    class(psb_i_base_vect_type), optional :: imold
    integer(psb_ipk_), optional :: partition, nrl,iv(:)

    ! Local variables.

    integer(psb_ipk_), parameter :: nb=20
    type(psb_d_csc_sparse_mat)  :: acsc
    type(psb_d_coo_sparse_mat)  :: acoo
    type(psb_d_csr_sparse_mat)  :: acsr
    real(psb_dpk_)           :: zt(nb),x,y,z
    integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_
    integer(psb_lpk_) :: m,n,glob_row,nt
    integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
    ! For 3D partition
    ! Note: integer control variables going directly into an MPI call
    ! must be 4 bytes, i.e. psb_mpk_
    integer(psb_mpk_) :: npdims(3), npp, minfo
    integer(psb_ipk_) :: npx,npy,npz, iamx,iamy,iamz,mynx,myny,mynz
    integer(psb_ipk_), allocatable :: bndx(:),bndy(:),bndz(:)
    ! Process grid
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: icoeff
    integer(psb_lpk_), allocatable     :: irow(:),icol(:),myidx(:)
    real(psb_dpk_), allocatable :: val(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_dpk_)            :: deltah, sqdeltah, deltah2
    real(psb_dpk_), parameter :: rhs=dzero,one=done,zero=dzero
    real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    procedure(d_func_3d), pointer :: f_
    character(len=20)  :: name, ch_err,tmpfmt

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)


    if (present(f)) then 
      f_ => f
    else
      f_ => d_null_func_3d
    end if

    deltah   = done/(idim+1)
    sqdeltah = deltah*deltah
    deltah2  = (2*done)* deltah

    if (present(partition)) then
      if ((1<= partition).and.(partition <= 3)) then
        partition_ = partition
      else
        write(*,*) 'Invalid partition choice ',partition,' defaulting to 3'
        partition_ = 3
      end if
    else
      partition_ = 3
    end if
    
    ! initialize array descriptor and sparse matrix storage. provide an
    ! estimate of the number of non zeroes 
    
    m   = (1_psb_lpk_*idim)*idim*idim
    n   = m
    nnz = 7*((n+np-1)/np)
    if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n
    t0 = psb_wtime()
    select case(partition_)
    case(1)
      ! A BLOCK partition 
      if (present(nrl)) then 
        nr = nrl
      else
        !
        ! Using a simple BLOCK distribution.
        !
        nt = (m+np-1)/np
        nr = max(0,min(nt,m-(iam*nt)))
      end if

      nt = nr
      call psb_sum(ctxt,nt) 
      if (nt /= m) then 
        write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return    
      end if

      !
      ! First example  of use of CDALL: specify for each process a number of
      ! contiguous rows
      ! 
      call psb_cdall(ctxt,desc_a,info,nl=nr)
      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

    case(2)
      ! A  partition  defined by the user through IV
      
      if (present(iv)) then 
        if (size(iv) /= m) then
          write(psb_err_unit,*) iam, 'Initialization error: wrong IV size',size(iv),m
          info = -1
          call psb_barrier(ctxt)
          call psb_abort(ctxt)
          return    
        end if
      else
        write(psb_err_unit,*) iam, 'Initialization error: IV not present'
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return    
      end if

      !
      ! Second example  of use of CDALL: specify for each row the
      ! process that owns it 
      ! 
      call psb_cdall(ctxt,desc_a,info,vg=iv)
      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

    case(3)
      ! A 3-dimensional partition

      ! A nifty MPI function will split the process list
      npdims = 0
      call mpi_dims_create(np,3,npdims,info)
      npx = npdims(1)
      npy = npdims(2)
      npz = npdims(3)

      allocate(bndx(0:npx),bndy(0:npy),bndz(0:npz))
      ! We can reuse idx2ijk for process indices as well. 
      call idx2ijk(iamx,iamy,iamz,iam,npx,npy,npz,base=0)
      ! Now let's split the 3D cube in hexahedra
      call dist1Didx(bndx,idim,npx)
      mynx = bndx(iamx+1)-bndx(iamx)
      call dist1Didx(bndy,idim,npy)
      myny = bndy(iamy+1)-bndy(iamy)
      call dist1Didx(bndz,idim,npz)
      mynz = bndz(iamz+1)-bndz(iamz)

      ! How many indices do I own? 
      nlr = mynx*myny*mynz
      allocate(myidx(nlr))
      ! Now, let's generate the list of indices I own
      nr = 0
      do i=bndx(iamx),bndx(iamx+1)-1
        do j=bndy(iamy),bndy(iamy+1)-1
          do k=bndz(iamz),bndz(iamz+1)-1
            nr = nr + 1
            call ijk2idx(myidx(nr),i,j,k,idim,idim,idim)
          end do
        end do
      end do
      if (nr /= nlr) then
        write(psb_err_unit,*) iam,iamx,iamy,iamz, 'Initialization error: NR vs NLR ',&
             & nr,nlr,mynx,myny,mynz
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
      end if

      !
      ! Third example  of use of CDALL: specify for each process
      ! the set of global indices it owns.
      ! 
      call psb_cdall(ctxt,desc_a,info,vl=myidx)
      
    case default
      write(psb_err_unit,*) iam, 'Initialization error: should not get here'
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end select

    
    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
    ! define  rhs from boundary conditions; also build initial guess 
    if (info == psb_success_) call psb_geall(xv,desc_a,info)
    if (info == psb_success_) call psb_geall(bv,desc_a,info)

    call psb_barrier(ctxt)
    talc = psb_wtime()-t0

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    ! we build an auxiliary matrix consisting of one row at a
    ! time; just a small matrix. might be extended to generate 
    ! a bunch of rows per call. 
    ! 
    allocate(val(20*nb),irow(20*nb),&
         &icol(20*nb),stat=info)
    if (info /= psb_success_ ) then 
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif


    ! loop over rows belonging to current process in a block
    ! distribution.

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    do ii=1, nlr,nb
      ib = min(nb,nlr-ii+1) 
      icoeff = 1
      do k=1,ib
        i=ii+k-1
        ! local matrix pointer 
        glob_row=myidx(i)
        ! compute gridpoint coordinates
        call idx2ijk(ix,iy,iz,glob_row,idim,idim,idim)
        ! x, y, z coordinates
        x = (ix-1)*deltah
        y = (iy-1)*deltah
        z = (iz-1)*deltah
        zt(k) = f_(x,y,z)
        ! internal point: build discretization
        !   
        !  term depending on   (x-1,y,z)
        !
        val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
        if (ix == 1) then 
          zt(k) = g(dzero,y,z)*(-val(icoeff)) + zt(k)
        else
          call ijk2idx(icol(icoeff),ix-1,iy,iz,idim,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y-1,z)
        val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
        if (iy == 1) then 
          zt(k) = g(x,dzero,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy-1,iz,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y,z-1)
        val(icoeff)=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
        if (iz == 1) then 
          zt(k) = g(x,y,dzero)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy,iz-1,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif

        !  term depending on     (x,y,z)
        val(icoeff)=(2*done)*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah &
             & + c(x,y,z)
        call ijk2idx(icol(icoeff),ix,iy,iz,idim,idim,idim)          
        irow(icoeff) = glob_row
        icoeff       = icoeff+1                  
        !  term depending on     (x,y,z+1)
        val(icoeff)=-a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
        if (iz == idim) then 
          zt(k) = g(x,y,done)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy,iz+1,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y+1,z)
        val(icoeff)=-a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
        if (iy == idim) then 
          zt(k) = g(x,done,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy+1,iz,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x+1,y,z)
        val(icoeff)=-a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
        if (ix==idim) then 
          zt(k) = g(done,y,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix+1,iy,iz,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif

      end do
      call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
      if(info /= psb_success_) exit
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
      if(info /= psb_success_) exit
      zt(:)=dzero
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
      if(info /= psb_success_) exit
    end do

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(val,irow,icol)

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info,mold=imold)
    tcdasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    if (info == psb_success_) then 
      if (present(amold)) then 
        call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,mold=amold)
      else
        call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      end if
    end if
    call psb_barrier(ctxt)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (info == psb_success_) call psb_geasb(xv,desc_a,info,mold=vmold)
    if (info == psb_success_) call psb_geasb(bv,desc_a,info,mold=vmold)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    tasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    ttot = psb_wtime() - t0 

    call psb_amx(ctxt,talc)
    call psb_amx(ctxt,tgen)
    call psb_amx(ctxt,tasb)
    call psb_amx(ctxt,ttot)
    if(iam == psb_root_) then
      tmpfmt = a%get_fmt()
      write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")')&
           &   tmpfmt
      write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
      write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
      write(psb_out_unit,'("-desc asbly  time : ",es12.5)') tcdasb
      write(psb_out_unit,'("- mat asbly  time : ",es12.5)') tasb
      write(psb_out_unit,'("-total       time : ",es12.5)') ttot

    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine amg_d_gen_pde3d

end module amg_d_pde3d_mod

program amg_d_pde3d
  use psb_base_mod
  use amg_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use data_input
  use amg_d_pde3d_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim
  integer(psb_epk_) :: system_size

  ! miscellaneous 
  real(psb_dpk_) :: t1, t2, tprec, thier, tslv

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(amg_dprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: x,b,r
  ! parallel environment
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, resmx, resmxp

  ! Krylov solver data
  type solverdata
    character(len=40)  :: kmethd      ! Krylov solver
    integer(psb_ipk_)  :: istopc      ! stopping criterion
    integer(psb_ipk_)  :: itmax       ! maximum number of iterations
    integer(psb_ipk_)  :: itrace      ! tracing
    integer(psb_ipk_)  :: irst        ! restart
    real(psb_dpk_)     :: eps         ! stopping tolerance
  end type solverdata
  type(solverdata)       :: s_choice

  ! preconditioner data
  type precdata

    ! preconditioner type
    character(len=40)  :: descr       ! verbose description of the prec
    character(len=10)  :: ptype       ! preconditioner type

    integer(psb_ipk_)  :: outer_sweeps ! number of outer sweeps: sweeps for 1-level,
                                       ! AMG cycles for ML
    ! general AMG data
    character(len=16)  :: mlcycle      ! AMG cycle type
    integer(psb_ipk_)  :: maxlevs     ! maximum number of levels in AMG preconditioner

    ! AMG aggregation
    character(len=16)  :: aggr_prol    ! aggregation type: SMOOTHED, NONSMOOTHED
    character(len=16)  :: par_aggr_alg    ! parallel aggregation algorithm: DEC, SYMDEC
    character(len=16)  :: aggr_ord    ! ordering for aggregation: NATURAL, DEGREE
    character(len=16)  :: aggr_filter ! filtering: FILTER, NO_FILTER
    real(psb_dpk_)     :: mncrratio  ! minimum aggregation ratio
    real(psb_dpk_), allocatable :: athresv(:) ! smoothed aggregation threshold vector
    integer(psb_ipk_)  :: thrvsz      ! size of threshold vector
    real(psb_dpk_)     :: athres      ! smoothed aggregation threshold
    integer(psb_ipk_)  :: csize       ! minimum size of coarsest matrix

    ! AMG smoother or pre-smoother; also 1-lev preconditioner
    character(len=16)  :: smther      ! (pre-)smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps     ! (pre-)smoother / 1-lev prec. sweeps
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    character(len=16)  :: restr       ! restriction over application of AS
    character(len=16)  :: prol        ! prolongation over application of AS
    character(len=16)  :: solve       ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: fill        ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: thr         ! threshold for ILUT factorization

    ! AMG post-smoother; ignored by 1-lev preconditioner
    character(len=16)  :: smther2     ! post-smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps2    ! post-smoother sweeps
    integer(psb_ipk_)  :: novr2       ! number of overlap layers
    character(len=16)  :: restr2      ! restriction  over application of AS
    character(len=16)  :: prol2       ! prolongation over application of AS
    character(len=16)  :: solve2      ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: fill2       ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: thr2        ! threshold for ILUT factorization

    ! coarsest-level solver
    character(len=16)  :: cmat        ! coarsest matrix layout: REPL, DIST
    character(len=16)  :: csolve      ! coarsest-lev solver: BJAC, SLUDIST (distr.
                                      ! mat.); UMF, MUMPS, SLU, ILU, ILUT, MILU
                                      ! (repl. mat.)
    character(len=16)  :: csbsolve    ! coarsest-lev local subsolver: ILU, ILUT,
                                      ! MILU, UMF, MUMPS, SLU
    integer(psb_ipk_)  :: cfill       ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: cthres      ! threshold for ILUT factorization
    integer(psb_ipk_)  :: cjswp       ! sweeps for GS or JAC coarsest-lev subsolver

  end type precdata
  type(precdata)       :: p_choice

  ! other variables
  integer(psb_ipk_)  :: info, i, k
  character(len=20)  :: name,ch_err

  info=psb_success_


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='amg_d_pde3d'
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to MLD2P4 version: ',amg_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  !
  !  get parameters
  !
  call get_parms(ctxt,afmt,idim,s_choice,p_choice)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,info)  
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='amg_gen_pde3d'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')
  !
  ! initialize the preconditioner
  !
  call prec%init(ctxt,p_choice%ptype,info)
  select case(trim(psb_toupper(p_choice%ptype)))
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','L1-JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    
  case ('BJAC')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case('AS')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_ovr',         p_choice%novr,    info)
    call prec%set('sub_restr',       p_choice%restr,   info)
    call prec%set('sub_prol',        p_choice%prol,    info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)
    
  case ('ML') 
    ! multilevel preconditioner

    call prec%set('ml_cycle',        p_choice%mlcycle,    info)
    call prec%set('outer_sweeps',    p_choice%outer_sweeps,info)
    if (p_choice%csize>0)&
         & call prec%set('min_coarse_size', p_choice%csize,      info)
    if (p_choice%mncrratio>1)&
         & call prec%set('min_cr_ratio',   p_choice%mncrratio, info)
    if (p_choice%maxlevs>0)&
         & call prec%set('max_levs',    p_choice%maxlevs,    info)
    if (p_choice%athres >= dzero) &
         & call prec%set('aggr_thresh',     p_choice%athres,  info)
    if (p_choice%thrvsz>0) then
      do k=1,min(p_choice%thrvsz,size(prec%precv)-1)
        call prec%set('aggr_thresh',     p_choice%athresv(k),  info,ilev=(k+1))
      end do
    end if

    call prec%set('aggr_prol',       p_choice%aggr_prol,   info)
    call prec%set('par_aggr_alg',    p_choice%par_aggr_alg,   info)
    call prec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call prec%set('aggr_filter',     p_choice%aggr_filter,info)


    call prec%set('smoother_type',   p_choice%smther,     info)
    call prec%set('smoother_sweeps', p_choice%jsweeps,    info)

    select case (psb_toupper(p_choice%smther))
    case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI')
      ! do nothing
    case default
      call prec%set('sub_ovr',         p_choice%novr,       info)
      call prec%set('sub_restr',       p_choice%restr,      info)
      call prec%set('sub_prol',        p_choice%prol,       info)
      call prec%set('sub_solve',       p_choice%solve,      info)
      call prec%set('sub_fillin',      p_choice%fill,       info)
      call prec%set('sub_iluthrs',     p_choice%thr,        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call prec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call prec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
      select case (psb_toupper(p_choice%smther2))
      case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI')
        ! do nothing
      case default
        call prec%set('sub_ovr',         p_choice%novr2,     info,pos='post')
        call prec%set('sub_restr',       p_choice%restr2,    info,pos='post')
        call prec%set('sub_prol',        p_choice%prol2,     info,pos='post')
        call prec%set('sub_solve',       p_choice%solve2,    info,pos='post')
        call prec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call prec%set('sub_iluthrs',     p_choice%thr2,      info,pos='post')
      end select
    end if

    call prec%set('coarse_solve',    p_choice%csolve,    info)
    if (psb_toupper(p_choice%csolve) == 'BJAC') &
         &  call prec%set('coarse_subsolve', p_choice%csbsolve,  info)
    call prec%set('coarse_mat',      p_choice%cmat,      info)
    call prec%set('coarse_fillin',   p_choice%cfill,     info)
    call prec%set('coarse_iluthrs',  p_choice%cthres,    info)
    call prec%set('coarse_sweeps',   p_choice%cjswp,     info)

  end select
  
  ! build the preconditioner
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  thier = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_hierarchy_bld')
    goto 9999
  end if
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info)
  tprec = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ctxt, thier)
  call psb_amx(ctxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')thier+tprec
    write(psb_out_unit,'(" ")')
  end if

  !
  ! iterative method parameters 
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_krylov(s_choice%kmethd,a,prec,b,x,s_choice%eps,&
       & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,itrace=s_choice%itrace,&
       & istop=s_choice%istopc,irst=s_choice%irst)
  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1

  call psb_amx(ctxt,tslv)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1
  call psb_amx(ctxt,tslv)

  ! compute residual norms
  call psb_geall(r,desc_a,info)
  call r%zero()
  call psb_geasb(r,desc_a,info)
  call psb_geaxpby(done,b,dzero,r,desc_a,info)
  call psb_spmm(-done,a,x,done,r,desc_a,info)
  resmx  = psb_genrm2(r,desc_a,info)
  resmxp = psb_geamax(r,desc_a,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)
  call prec%descr(iout=psb_out_unit)
  if (iam == psb_root_) then 
    write(psb_out_unit,'("Computed solution on ",i8," processors")')  np
    write(psb_out_unit,'("Linear system size                 : ",i12)') system_size
    write(psb_out_unit,'("Krylov method                      : ",a)') trim(s_choice%kmethd)
    write(psb_out_unit,'("Preconditioner                     : ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Iterations to convergence          : ",i12)')    iter
    write(psb_out_unit,'("Relative error estimate on exit    : ",es12.5)') err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)')    prec%get_nlevs()
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)') thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)') tprec
    write(psb_out_unit,'("Total time for preconditioner      : ",es12.5)') tprec+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)') tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)') tslv/iter
    write(psb_out_unit,'("Total time                         : ",es12.5)') tslv+tprec+thier
    write(psb_out_unit,'("Residual 2-norm                    : ",es12.5)') resmx
    write(psb_out_unit,'("Residual inf-norm                  : ",es12.5)') resmxp
    write(psb_out_unit,'("Total memory occupation for A      : ",i12)') amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i12)') descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i12)') precsize
    write(psb_out_unit,'("Storage format for A               : ",a  )') a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )') desc_a%get_fmt()

  end if

  !  
  !  cleanup storage and exit
  !
  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get iteration parameters from standard input
  !
  !
  ! get iteration parameters from standard input
  !
  subroutine get_parms(ctxt,afmt,idim,solve,prec)

    implicit none

    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: idim
    character(len=*)    :: afmt
    type(solverdata)    :: solve
    type(precdata)      :: prec
    integer(psb_ipk_)   :: iam, nm, np, inp_unit
    character(len=1024)   :: filename

    call psb_info(ctxt,iam,np)

    if (iam == psb_root_) then
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ctxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      ! read input data
      !
      call read_data(afmt,inp_unit)            ! matrix storage format
      call read_data(idim,inp_unit)            ! Discretization grid size
      ! Krylov solver data
      call read_data(solve%kmethd,inp_unit)    ! Krylov solver
      call read_data(solve%istopc,inp_unit)    ! stopping criterion
      call read_data(solve%itmax,inp_unit)     ! max num iterations
      call read_data(solve%itrace,inp_unit)    ! tracing
      call read_data(solve%irst,inp_unit)      ! restart
      call read_data(solve%eps,inp_unit)       ! tolerance
      ! preconditioner type
      call read_data(prec%descr,inp_unit)      ! verbose description of the prec
      call read_data(prec%ptype,inp_unit)      ! preconditioner type
      ! First smoother / 1-lev preconditioner
      call read_data(prec%smther,inp_unit)     ! smoother type
      call read_data(prec%jsweeps,inp_unit)    ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%novr,inp_unit)       ! number of overlap layers
      call read_data(prec%restr,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve,inp_unit)      ! local subsolver
      call read_data(prec%fill,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%thr,inp_unit)        ! threshold for ILUT
      ! Second smoother/ AMG post-smoother (if NONE ignored in main)
      call read_data(prec%smther2,inp_unit)     ! smoother type
      call read_data(prec%jsweeps2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%novr2,inp_unit)       ! number of overlap layers
      call read_data(prec%restr2,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol2,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve2,inp_unit)      ! local subsolver
      call read_data(prec%fill2,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%thr2,inp_unit)        ! threshold for ILUT
      ! general AMG data
      call read_data(prec%mlcycle,inp_unit)     ! AMG cycle type
      call read_data(prec%outer_sweeps,inp_unit) ! number of 1lev/outer sweeps
      call read_data(prec%maxlevs,inp_unit)    ! max number of levels in AMG prec
      call read_data(prec%csize,inp_unit)       ! min size coarsest mat
      ! aggregation
      call read_data(prec%aggr_prol,inp_unit)    ! aggregation type
      call read_data(prec%par_aggr_alg,inp_unit)    ! parallel aggregation alg
      call read_data(prec%aggr_ord,inp_unit)    ! ordering for aggregation
      call read_data(prec%aggr_filter,inp_unit) ! filtering
      call read_data(prec%mncrratio,inp_unit)  ! minimum aggregation ratio
      call read_data(prec%thrvsz,inp_unit)      ! size of aggr thresh vector
      if (prec%thrvsz > 0) then
        call psb_realloc(prec%thrvsz,prec%athresv,info)
        call read_data(prec%athresv,inp_unit)   ! aggr thresh vector
      else
        read(inp_unit,*)                        ! dummy read to skip a record
      end if
      call read_data(prec%athres,inp_unit)      ! smoothed aggr thresh
      ! coasest-level solver
      call read_data(prec%csolve,inp_unit)      ! coarsest-lev solver
      call read_data(prec%csbsolve,inp_unit)    ! coarsest-lev subsolver
      call read_data(prec%cmat,inp_unit)        ! coarsest mat layout
      call read_data(prec%cfill,inp_unit)       ! fill-in for incompl LU
      call read_data(prec%cthres,inp_unit)      ! Threshold for ILUT
      call read_data(prec%cjswp,inp_unit)       ! sweeps for GS/JAC subsolver
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    call psb_bcast(ctxt,afmt)
    call psb_bcast(ctxt,idim)

    call psb_bcast(ctxt,solve%kmethd)
    call psb_bcast(ctxt,solve%istopc)
    call psb_bcast(ctxt,solve%itmax)
    call psb_bcast(ctxt,solve%itrace)
    call psb_bcast(ctxt,solve%irst)
    call psb_bcast(ctxt,solve%eps)

    call psb_bcast(ctxt,prec%descr)
    call psb_bcast(ctxt,prec%ptype)

    ! broadcast first (pre-)smoother / 1-lev prec data
    call psb_bcast(ctxt,prec%smther)     
    call psb_bcast(ctxt,prec%jsweeps)
    call psb_bcast(ctxt,prec%novr)
    call psb_bcast(ctxt,prec%restr)
    call psb_bcast(ctxt,prec%prol)
    call psb_bcast(ctxt,prec%solve)
    call psb_bcast(ctxt,prec%fill)
    call psb_bcast(ctxt,prec%thr)
    ! broadcast second (post-)smoother 
    call psb_bcast(ctxt,prec%smther2)
    call psb_bcast(ctxt,prec%jsweeps2)
    call psb_bcast(ctxt,prec%novr2)
    call psb_bcast(ctxt,prec%restr2)
    call psb_bcast(ctxt,prec%prol2)
    call psb_bcast(ctxt,prec%solve2)
    call psb_bcast(ctxt,prec%fill2)
    call psb_bcast(ctxt,prec%thr2)
    
    ! broadcast AMG parameters
    call psb_bcast(ctxt,prec%mlcycle)
    call psb_bcast(ctxt,prec%outer_sweeps)
    call psb_bcast(ctxt,prec%maxlevs)
    
    call psb_bcast(ctxt,prec%aggr_prol)
    call psb_bcast(ctxt,prec%par_aggr_alg)
    call psb_bcast(ctxt,prec%aggr_ord)
    call psb_bcast(ctxt,prec%aggr_filter)
    call psb_bcast(ctxt,prec%mncrratio)
    call psb_bcast(ctxt,prec%thrvsz)
    if (prec%thrvsz > 0) then
      if (iam /= psb_root_) call psb_realloc(prec%thrvsz,prec%athresv,info)
      call psb_bcast(ctxt,prec%athresv)
    end if
    call psb_bcast(ctxt,prec%athres)
    
    call psb_bcast(ctxt,prec%csize)
    call psb_bcast(ctxt,prec%cmat)
    call psb_bcast(ctxt,prec%csolve)
    call psb_bcast(ctxt,prec%csbsolve)
    call psb_bcast(ctxt,prec%cfill)
    call psb_bcast(ctxt,prec%cthres)
    call psb_bcast(ctxt,prec%cjswp)


  end subroutine get_parms

end program amg_d_pde3d

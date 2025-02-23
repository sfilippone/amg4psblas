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
! File: amg_dexample_cuda.f90
!
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. The solver is CG, coupled with one of the
! following multi-level preconditioner, as explained in Section 4.2 of
! the AMG4PSBLAS User's and Reference Guide:
!
! - choice = 1, a V-cycle with decoupled smoothed aggregation, 4 Jacobi
! sweeps as pre/post-smoother and 8 Jacobi sweeps as coarsest-level
! solver with replicated coarsest matrix
!
! - choice = 2, a W-cycle based on the coupled aggregation relying on matching,
! with maximum size of aggregates equal to 8 and smoothed prolongators,
! 2 sweeps of Block-Jacobi ipre/post-smoother using approximate inverse INVK and
! 4 sweeps of Block-Jacobi with INVK as coarsest-level solver on distributed
! coarsest matrix 
!
! The matrix and the rhs are read from files (if an rhs is not available, the
! unit rhs is set).
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
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
module amg_d_pde3d_poisson_mod
  use psb_base_mod, only : psb_dpk_, done, dzero
  real(psb_dpk_), save, private :: epsilon=done/80
contains
  subroutine pde_set_parm3d_poisson(dat)
    real(psb_dpk_), intent(in) :: dat
    epsilon = dat
  end subroutine pde_set_parm3d_poisson
  !
  ! functions parametrizing the differential equation
  !
  function b1_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) :: b1_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    b1_poisson=dzero
  end function b1_poisson
  function b2_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  b2_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    b2_poisson=dzero
  end function b2_poisson
  function b3_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  b3_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    b3_poisson=dzero
  end function b3_poisson
  function c_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  c_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    c_poisson=dzero
  end function c_poisson
  function a1_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a1_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    a1_poisson=epsilon
  end function a1_poisson
  function a2_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a2_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    a2_poisson=epsilon
  end function a2_poisson
  function a3_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  a3_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    a3_poisson=epsilon
  end function a3_poisson
  function g_poisson(x,y,z)
    implicit none 
    real(psb_dpk_) ::  g_poisson
    real(psb_dpk_), intent(in) :: x,y,z
    g_poisson = dzero
    if (x == done) then
      g_poisson = done
    else if (x == dzero) then
      g_poisson = done
    end if
  end function g_poisson
end module amg_d_pde3d_poisson_mod

module amg_d_genpde_mod

  use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_desc_type,&
       &  psb_dspmat_type, psb_d_vect_type, dzero, done,&
       &  psb_d_base_sparse_mat, psb_d_base_vect_type, psb_i_base_vect_type

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

  interface
    function d_func_2d(x,y) result(val)
      import :: psb_dpk_
      real(psb_dpk_), intent(in) :: x,y
      real(psb_dpk_) :: val
    end function d_func_2d
  end interface

  interface amg_gen_pde2d
    module procedure  amg_d_gen_pde2d
  end interface amg_gen_pde2d

contains

  function d_null_func_2d(x,y) result(val)

    real(psb_dpk_), intent(in) :: x,y
    real(psb_dpk_) :: val

    val = dzero

  end function d_null_func_2d

  function d_null_func_3d(x,y,z) result(val)

    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val

    val = dzero

  end function d_null_func_3d

  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs.
  !
  subroutine amg_d_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info,f,amold,vmold,partition, nrl,iv)
    use psb_base_mod
    use psb_util_mod
    !
    !   Discretizes the partial differential equation
    !
    !    d a1 d(u)  d a1 d(u)   d a1 d(u)    b1 d(u)   b2 d(u)    b3 d(u)
    ! -   ------ -  ------   -  ------     +  -----  +  ------  +  ------ + c u = f
    !     dx  dx     dy  dy      dz  dz         dx        dy         dz
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
    procedure(d_func_3d)  :: b1,b2,b3,c,a1,a2,a3,g
    integer(psb_ipk_)     :: idim
    type(psb_dspmat_type) :: a
    type(psb_d_vect_type) :: xv,bv
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_)     :: info
    type(psb_ctxt_type)   :: ctxt
    character             :: afmt*5
    procedure(d_func_3d), optional :: f
    class(psb_d_base_sparse_mat), optional :: amold
    class(psb_d_base_vect_type), optional :: vmold
    integer(psb_ipk_), optional :: partition, nrl,iv(:)

    ! Local variables.

    integer(psb_ipk_), parameter :: nb=20
    type(psb_d_csc_sparse_mat)  :: acsc
    type(psb_d_coo_sparse_mat)  :: acoo
    type(psb_d_csr_sparse_mat)  :: acsr
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
    integer(psb_lpk_), allocatable     :: myidx(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_dpk_)            :: deltah, sqdeltah, deltah2
    real(psb_dpk_), parameter :: rhs=dzero,one=done,zero=dzero
    real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    procedure(d_func_3d), pointer :: f_
    character(len=20)  :: name, ch_err,tmpfmt

    info = psb_success_
    name = 'd_create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)


    if (present(f)) then
      f_ => f
    else
      f_ => d_null_func_3d
    end if

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
    deltah   = done/(idim+2)
    sqdeltah = deltah*deltah
    deltah2  = 2.0_psb_dpk_* deltah

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

      !
      ! Specify process topology
      !
      block
        !
        ! Use adjcncy methods 
        ! 
        integer(psb_mpk_), allocatable :: neighbours(:)
        integer(psb_mpk_) :: cnt
        logical, parameter :: debug_adj=.true.
        if (debug_adj.and.(np > 1)) then 
          cnt = 0
          allocate(neighbours(np))
          if (iamx < npx-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx+1,iamy,iamz,npx,npy,npz,base=0)
          end if
          if (iamy < npy-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy+1,iamz,npx,npy,npz,base=0)
          end if
          if (iamz < npz-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy,iamz+1,npx,npy,npz,base=0)
          end if
          if (iamx >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx-1,iamy,iamz,npx,npy,npz,base=0)
          end if
          if (iamy >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy-1,iamz,npx,npy,npz,base=0)
          end if
          if (iamz >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy,iamz-1,npx,npy,npz,base=0)
          end if
          call psb_realloc(cnt, neighbours,info)
          call desc_a%set_p_adjcncy(neighbours)
          !write(0,*) iam,' Check on neighbours: ',desc_a%get_p_adjcncy()
        end if
      end block
      
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

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    !$omp parallel shared(deltah,myidx,a,desc_a)
    !
    block 
      integer(psb_ipk_) :: i,j,k,ii,ib,icoeff, ix,iy,iz, ith,nth
      integer(psb_lpk_) :: glob_row
      integer(psb_lpk_), allocatable :: irow(:),icol(:)
      real(psb_dpk_), allocatable :: val(:)
      real(psb_dpk_)     :: x,y,z, zt(nb)
      nth = 1
      ith = 0
      allocate(val(20*nb),irow(20*nb),&
           &icol(20*nb),stat=info)
      if (info /= psb_success_ ) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        !goto 9999
      endif

      !$omp  do schedule(dynamic)
      !     
      do ii=1, nlr, nb
        if (info /= psb_success_) cycle
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
        !write(0,*) ' Outer in_parallel ',omp_in_parallel()
        call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
        if(info /= psb_success_) cycle
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
        if(info /= psb_success_) cycle
        zt(:)=dzero
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
        if(info /= psb_success_) cycle
      end do
      !$omp end do

      deallocate(val,irow,icol)
    end block
    !$omp end parallel

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info)
    tcdasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    if (info == psb_success_) then
      if (present(amold)) then
        call psb_spasb(a,desc_a,info,mold=amold)
      else
        call psb_spasb(a,desc_a,info,afmt=afmt)
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



  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs.
  !
  subroutine amg_d_gen_pde2d(ctxt,idim,a,bv,xv,desc_a,afmt,&
       & a1,a2,b1,b2,c,g,info,f,amold,vmold,partition, nrl,iv)
    use psb_base_mod
    use psb_util_mod
    !
    !   Discretizes the partial differential equation
    !
    !    d    d(u)   d     d(u)    b1 d(u)  b2 d(u)
    ! - -- a1 ---- - -- a1 ---- +  -----  +  ------  + c u = f
    !    dx   dx     dy    dy       dx       dy
    !
    ! with Dirichlet boundary conditions
    !   u = g
    !
    !  on the unit square  0<=x,y<=1.
    !
    !
    ! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
    !
    implicit none
    procedure(d_func_2d)  :: b1,b2,c,a1,a2,g
    integer(psb_ipk_)     :: idim
    type(psb_dspmat_type) :: a
    type(psb_d_vect_type) :: xv,bv
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_)     :: info
    type(psb_ctxt_type)   :: ctxt
    character             :: afmt*5
    procedure(d_func_2d), optional :: f
    class(psb_d_base_sparse_mat), optional :: amold
    class(psb_d_base_vect_type), optional :: vmold
    integer(psb_ipk_), optional :: partition, nrl,iv(:)
    ! Local variables.

    integer(psb_ipk_), parameter :: nb=20
    type(psb_d_csc_sparse_mat)  :: acsc
    type(psb_d_coo_sparse_mat)  :: acoo
    type(psb_d_csr_sparse_mat)  :: acsr
    integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_
    integer(psb_lpk_) :: m,n,glob_row,nt
    integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
    ! For 2D partition
    ! Note: integer control variables going directly into an MPI call
    ! must be 4 bytes, i.e. psb_mpk_
    integer(psb_mpk_) :: npdims(2), npp, minfo
    integer(psb_ipk_) :: npx,npy,iamx,iamy,mynx,myny
    integer(psb_ipk_), allocatable :: bndx(:),bndy(:)
    ! Process grid
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: icoeff
    integer(psb_lpk_), allocatable     :: myidx(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_dpk_)            :: deltah, sqdeltah, deltah2, dd
    real(psb_dpk_), parameter :: rhs=0.d0,one=done,zero=0.d0
    real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    procedure(d_func_2d), pointer :: f_
    character(len=20)  :: name, ch_err,tmpfmt

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)


    if (present(f)) then
      f_ => f
    else
      f_ => d_null_func_2d
    end if

    deltah   = done/(idim+2)
    sqdeltah = deltah*deltah
    deltah2  = 2.0_psb_dpk_* deltah


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

    m   = (1_psb_lpk_)*idim*idim
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
      ! A 2-dimensional partition

      ! A nifty MPI function will split the process list
      npdims = 0
      call mpi_dims_create(np,2,npdims,info)
      npx = npdims(1)
      npy = npdims(2)

      allocate(bndx(0:npx),bndy(0:npy))
      ! We can reuse idx2ijk for process indices as well.
      call idx2ijk(iamx,iamy,iam,npx,npy,base=0)
      ! Now let's split the 2D square in rectangles
      call dist1Didx(bndx,idim,npx)
      mynx = bndx(iamx+1)-bndx(iamx)
      call dist1Didx(bndy,idim,npy)
      myny = bndy(iamy+1)-bndy(iamy)

      ! How many indices do I own?
      nlr = mynx*myny
      allocate(myidx(nlr))
      ! Now, let's generate the list of indices I own
      nr = 0
      do i=bndx(iamx),bndx(iamx+1)-1
        do j=bndy(iamy),bndy(iamy+1)-1
          nr = nr + 1
          call ijk2idx(myidx(nr),i,j,idim,idim)
        end do
      end do
      if (nr /= nlr) then
        write(psb_err_unit,*) iam,iamx,iamy, 'Initialization error: NR vs NLR ',&
             & nr,nlr,mynx,myny
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
      end if

      !
      ! Third example  of use of CDALL: specify for each process
      ! the set of global indices it owns.
      !
      call psb_cdall(ctxt,desc_a,info,vl=myidx)

      !
      ! Specify process topology
      !
      block
        !
        ! Use adjcncy methods 
        ! 
        integer(psb_mpk_), allocatable :: neighbours(:)
        integer(psb_mpk_) :: cnt
        logical, parameter :: debug_adj=.true.
        if (debug_adj.and.(np > 1)) then 
          cnt = 0
          allocate(neighbours(np))
          if (iamx < npx-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx+1,iamy,npx,npy,base=0)
          end if
          if (iamy < npy-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy+1,npx,npy,base=0)
          end if
          if (iamx >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx-1,iamy,npx,npy,base=0)
          end if
          if (iamy >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy-1,npx,npy,base=0)
          end if
          call psb_realloc(cnt, neighbours,info)
          call desc_a%set_p_adjcncy(neighbours)
          !write(0,*) iam,' Check on neighbours: ',desc_a%get_p_adjcncy()
        end if
      end block

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

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    !$omp parallel shared(deltah,myidx,a,desc_a)
    !
    block 
      integer(psb_ipk_) :: i,j,k,ii,ib,icoeff, ix,iy,iz, ith,nth
      integer(psb_lpk_) :: glob_row
      integer(psb_lpk_), allocatable :: irow(:),icol(:)
      real(psb_dpk_), allocatable :: val(:)
      real(psb_dpk_)     :: x,y,z, zt(nb)
      nth = 1
      ith = 0
      allocate(val(20*nb),irow(20*nb),&
           &icol(20*nb),stat=info)
      if (info /= psb_success_ ) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        !goto 9999
      endif

      ! loop over rows belonging to current process in a block
      ! distribution.
      !$omp  do schedule(dynamic)
      !     
      do ii=1, nlr,nb
        ib = min(nb,nlr-ii+1)
        icoeff = 1
        do k=1,ib
          i=ii+k-1
          ! local matrix pointer
          glob_row=myidx(i)
          ! compute gridpoint coordinates
          call idx2ijk(ix,iy,glob_row,idim,idim)
          ! x, y coordinates
          x = (ix-1)*deltah
          y = (iy-1)*deltah

          zt(k) = f_(x,y)
          ! internal point: build discretization
          !
          !  term depending on   (x-1,y)
          !
          val(icoeff) = -a1(x,y)/sqdeltah-b1(x,y)/deltah2
          if (ix == 1) then
            zt(k) = g(dzero,y)*(-val(icoeff)) + zt(k)
          else
            call ijk2idx(icol(icoeff),ix-1,iy,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y-1)
          val(icoeff)  = -a2(x,y)/sqdeltah-b2(x,y)/deltah2
          if (iy == 1) then
            zt(k) = g(x,dzero)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy-1,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif

          !  term depending on     (x,y)
          val(icoeff)=(2*done)*(a1(x,y) + a2(x,y))/sqdeltah + c(x,y)
          call ijk2idx(icol(icoeff),ix,iy,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
          !  term depending on     (x,y+1)
          val(icoeff)=-a2(x,y)/sqdeltah+b2(x,y)/deltah2
          if (iy == idim) then
            zt(k) = g(x,done)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy+1,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x+1,y)
          val(icoeff)=-a1(x,y)/sqdeltah+b1(x,y)/deltah2
          if (ix==idim) then
            zt(k) = g(done,y)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix+1,iy,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif

        end do
        call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
        if(info /= psb_success_) cycle
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
        if(info /= psb_success_) cycle
        zt(:)=dzero
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
        if(info /= psb_success_) cycle
      end do
      !$omp end do

      deallocate(val,irow,icol)
    end block
    !$omp end parallel

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info)
    tcdasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    if (info == psb_success_) then
      if (present(amold)) then
        call psb_spasb(a,desc_a,info,mold=amold)
      else
        call psb_spasb(a,desc_a,info,afmt=afmt)
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

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ctxt)
      return
    end if
    return
  end subroutine amg_d_gen_pde2d
end module amg_d_genpde_mod

program amg_dexample_cuda
  use psb_base_mod
  use amg_prec_mod
  use psb_linsolve_mod
  use psb_util_mod
#if defined(HAVE_CUDA)
  use psb_cuda_mod
#endif
  use data_input
  use amg_d_genpde_mod
  use amg_d_pde3d_poisson_mod
  implicit none

  ! input parameters

  ! sparse matrices
  type(psb_dspmat_type) :: A

  ! sparse matrices descriptor
  type(psb_desc_type):: desc_A

  ! preconditioner
  type(amg_dprec_type)  :: prec

  ! right-hand side, solution and residual vectors
  type(psb_d_vect_type) :: x, b, r
  ! GPU variables
  type(psb_d_csr_sparse_mat), target   :: acsr
  type(psb_d_base_vect_type), target   :: dvect
  type(psb_i_base_vect_type), target   :: ivect
#if defined(HAVE_CUDA)
  type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
  type(psb_d_cuda_hdiag_sparse_mat), target :: ahdiag
  type(psb_d_vect_cuda), target         :: dvgpu
  type(psb_i_vect_cuda), target         :: ivgpu
  class(psb_d_base_sparse_mat), pointer :: amold => ahlg
  class(psb_d_base_vect_type), pointer :: vmold => dvgpu
  class(psb_i_base_vect_type), pointer :: imold => ivgpu
#else
  class(psb_d_base_sparse_mat), pointer :: amold => acsr
  class(psb_d_base_vect_type), pointer :: vmold => dvect
  class(psb_i_base_vect_type), pointer :: imold => ivect
#endif
  ! solver and preconditioner parameters
  real(psb_dpk_)   :: tol, err
  integer          :: itmax, iter, istop
  integer          :: nlev

  ! parallel environment parameters
  type(psb_ctxt_type) :: ctxt
  integer             :: iam, np

  ! other variables
  integer            :: choice       
  integer            :: i,info,j
  integer(psb_epk_) :: amatsize, precsize, descsize
  integer(psb_epk_) :: system_size
  integer            :: idim, ierr, ircode
  real(psb_dpk_)     :: resmx, resmxp
  real(psb_dpk_)     :: t1, t2, tprec
  character(len=5)   :: afmt='CSR'
  character(len=20)  :: name, kmethod

  ! initialize the parallel environment

  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
#if defined(HAVE_CUDA)
  !
  ! BEWARE: if you have NGPUS  per node, the default is to
  ! attach to mod(IAM,NGPUS)
  !
  call psb_cuda_init(ctxt)
#endif

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif

  name='amg_dexample_cuda'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to AMG4PSBLAS version: ',amg_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
#if defined(HAVE_CUDA)
  write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
  write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_cuda_DeviceName())  
#endif
  ! get parameters

  call get_parms(ctxt,choice,idim,itmax,tol)

  !  allocate and fill in the coefficient matrix, rhs and initial guess

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
       & a1_poisson,a2_poisson,a3_poisson,&
       & b1_poisson,b2_poisson,b3_poisson,c_poisson,g_poisson,info)
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (iam == psb_root_) write(*,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(*,'(" ")')

  select case(choice)

 case(1)

    ! initialize a V-cycle preconditioner, relying on decoupled smoothed aggregation 
    ! with 4 Jacobi sweeps as pre/post-smoother
    ! and 8 Jacobi sweeps as coarsest-level solver on replicated coarsest matrix

    call prec%init(ctxt,'ML',info)
    call prec%set('SMOOTHER_TYPE','JACOBI',info)
    call prec%set('SMOOTHER_SWEEPS',4,info)
    call prec%set('COARSE_SOLVE','JACOBI',info)
    call prec%set('COARSE_SWEEPS',8,info)
    kmethod = 'CG'

  case(2)

   ! initialize a W-cycle preconditioner based on the coupled aggregation relying on matching,
   ! with maximum size of aggregates equal to 8 and smoothed prolongators,
   ! 2 sweeps of Block-Jacobi pre/post-smoother using approximate inverse INVK and
   ! 4 sweeps of Block-Jacobi with INVK on the coarsest level distributed matrix

    call prec%init(ctxt,'ML',info)
    call prec%set('PAR_AGGR_ALG','COUPLED',info)
    call prec%set('AGGR_TYPE','MATCHBOXP',info)
    call prec%set('AGGR_SIZE',8,info)
    call prec%set('ML_CYCLE','WCYCLE',info)
    call prec%set('SMOOTHER_TYPE','BJAC',info)
    call prec%set('SMOOTHER_SWEEPS',2,info)
    call prec%set('SUB_SOLVE','INVK',info)  
    call prec%set('COARSE_SOLVE','BJAC',info)
    call prec%set('COARSE_SUBSOLVE','INVK',info)
    call prec%set('COARSE_SWEEPS',4,info)
    call prec%set('COARSE_MAT','DIST',info)
    kmethod = 'CG'

  end select

  call psb_barrier(ctxt)
  t1 = psb_wtime()

  ! build the preconditioner
  call prec%hierarchy_build(A,desc_A,info)
  call prec%smoothers_build(A,desc_A,info, amold=amold, vmold=vmold, imold=imold)

  tprec = psb_wtime()-t1
  call psb_amx(ctxt, tprec)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_precbld')
    goto 9999
  end if

  ! set the solver parameters and the initial guess

  call psb_geall(x,desc_A,info)
  call x%zero()
  call psb_geasb(x,desc_A,info)

  ! Convert A, DESC_A,X,B to a GPU-enabled format
  call desc_a%cnv(mold=imold)
  call a%cscnv(info,mold=amold)
  call psb_geasb(x,desc_a,info,mold=vmold)
  call psb_geasb(b,desc_a,info,mold=vmold)


  ! solve Ax=b with preconditioned Krylov method

  call psb_barrier(ctxt)
  !
  ! Most preconditioners require auxiliary storage. When running
  ! on the HOST side, allocation/deallocation are usually very cheap
  ! and can be performed for every invocation of prec%apply.
  ! However when running on the DEVICE side, such memory management
  ! operations are global synchronization points, hence very costly.
  ! Thus the two methods below that preallocate the memory space
  ! prior to the invocation of the Krylov method, and release memory
  ! after the method has completed.
  !
  call prec%allocate_wrk(info,vmold=vmold)
  t1 = psb_wtime()

  call psb_krylov(kmethod,A,prec,b,x,tol,desc_A,info,&
       & itmax=itmax,iter=iter,err=err,itrace=1,istop=2)

  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)
  call prec%deallocate_wrk(info)

  call psb_geall(r,desc_A,info)
  call r%zero()
  call psb_geasb(r,desc_A,info)
  call psb_geaxpby(done,b,dzero,r,desc_A,info)
  call psb_spmm(-done,A,x,done,r,desc_A,info)
  resmx  = psb_genrm2(r,desc_A,info)
  resmxp = psb_geamax(r,desc_A,info)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)

  call prec%descr(info)

  if (iam == psb_root_) then
    write(*,'(" ")')
    write(*,'("Matrix from PDE example")')
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Linear system size        : ",i12)') system_size
    write(*,'("Krylov method             : ",a)') kmethod
    write(*,'("Iterations to convergence : ",i6)')iter
    write(*,'("Error estimate on exit    : ",es12.5)')err
    write(*,'("Time to build prec.       : ",es12.5)')tprec
    write(*,'("Time to solve system      : ",es12.5)')t2
    write(*,'("Time per iteration        : ",es12.5)')t2/(iter)
    write(*,'("Total time                : ",es12.5)')t2+tprec
    write(*,'("Residual 2-norm           : ",es12.5)')resmx
    write(*,'("Residual inf-norm         : ",es12.5)')resmxp
    write(*,'("Total memory occupation for A      : ",i12)')amatsize
    write(*,'("Total memory occupation for DESC_A : ",i12)')descsize
    write(*,'("Total memory occupation for PREC   : ",i12)')precsize
    write(*,'("Storage format for                A: ",a)')   trim(a%get_fmt())
    write(*,'("Storage format for           DESC_A: ",a)')   trim(desc_a%get_fmt())
  end if

  call psb_gefree(b, desc_A,info)
  call psb_gefree(x, desc_A,info)
  call psb_spfree(A, desc_A,info)
  call prec%free(info)
  call psb_cdfree(desc_A,info)
#if defined(HAVVE_CUDA)
  call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get parameters from standard input
  !
  subroutine get_parms(ctxt,choice,idim,itmax,tol)

    implicit none

    type(psb_ctxt_type) :: ctxt
    integer             :: choice, idim, itmax
    real(psb_dpk_)      :: tol
    integer             :: iam, np, inp_unit
    character(len=1024) :: filename

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
      ! read input parameters
      call read_data(choice,inp_unit)
      call read_data(idim,inp_unit)
      call read_data(itmax,inp_unit)
      call read_data(tol,inp_unit)
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if   
    end if

    call psb_bcast(ctxt,choice)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,itmax)
    call psb_bcast(ctxt,tol)

  end subroutine get_parms

end program amg_dexample_cuda

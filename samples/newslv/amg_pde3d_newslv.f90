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
      call idx2ijk(iamx,iamy,iamz,iam,npx,npy,npz,base=mzero)
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
        integer(psb_ipk_), allocatable :: neighbours(:)
        integer(psb_mpk_) :: cnt
        logical, parameter :: debug_adj=.true.
        if (debug_adj.and.(np > 1)) then 
          cnt = 0
          allocate(neighbours(np))
          if (iamx < npx-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx+1,iamy,iamz,npx,npy,npz,base=mzero)
          end if
          if (iamy < npy-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy+1,iamz,npx,npy,npz,base=mzero)
          end if
          if (iamz < npz-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy,iamz+1,npx,npy,npz,base=mzero)
          end if
          if (iamx >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx-1,iamy,iamz,npx,npy,npz,base=mzero)
          end if
          if (iamy >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy-1,iamz,npx,npy,npz,base=mzero)
          end if
          if (iamz >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy,iamz-1,npx,npy,npz,base=mzero)
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
      call idx2ijk(iamx,iamy,iam,npx,npy,base=mzero)
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
        integer(psb_ipk_), allocatable :: neighbours(:)
        integer(psb_mpk_) :: cnt
        logical, parameter :: debug_adj=.true.
        if (debug_adj.and.(np > 1)) then 
          cnt = 0
          allocate(neighbours(np))
          if (iamx < npx-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx+1,iamy,npx,npy,base=mzero)
          end if
          if (iamy < npy-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy+1,npx,npy,base=mzero)
          end if
          if (iamx >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx-1,iamy,npx,npy,base=mzero)
          end if
          if (iamy >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy-1,npx,npy,base=mzero)
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

program amg_d_pde3d
  use psb_base_mod
  use amg_prec_mod
  use psb_linsolve_mod
  use psb_util_mod
  use data_input
  use amg_d_genpde_mod
  use amg_d_pde3d_poisson_mod
  use amg_d_tlu_solver
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim

  ! miscellaneous 
  real(psb_dpk_) :: t1, t2, tprec, thier, tslv

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(amg_dprec_type)  :: prec
  type(amg_d_tlu_solver_type) :: tlusv
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: x,b
  ! parallel environment
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: iam, np

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps

  ! other variables
  integer(psb_ipk_)  :: info, i
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
    write(*,*) 'Welcome to AMG4PSBLAS version: ',amg_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  !
  !  get parameters
  !
  call get_parms(ctxt,kmethd,afmt,idim,istopc,itmax,itrace,irst,eps)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call amg_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,&
       & a1_poisson,a2_poisson,a3_poisson,&
       & b1_poisson,b2_poisson,b3_poisson,c_poisson,g_poisson,info)
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')
  !
  !  prepare the preconditioner: an ML with defaults, but with TLU solver at
  !  intermediate levels. All other parameters are at default values. 
  !  
  call prec%init(ctxt,'ML',       info)
  
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='prec%hierarchy_bld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  thier = psb_wtime()-t1
  nlv = prec%get_nlevs()
  call prec%set(tlusv,   info,ilev=1_psb_ipk_,ilmax=ione*max(1,nlv-1))
  
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='prec%smoothers_build'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  tprec = psb_wtime()-t1

  call psb_amx(ctxt,thier)
  call psb_amx(ctxt,tprec)

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Preconditioner time : ",es12.5)') tprec+thier
  if (iam == psb_root_) call prec%descr(info)
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')

  !
  ! iterative method parameters 
  !
  if(iam == psb_root_) &
       & write(psb_out_unit,'("Calling iterative method ",a)')kmethd
  call psb_barrier(ctxt)
  t1 = psb_wtime()  
  call psb_krylov(kmethd,a,prec,b,x,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)     

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1
  call psb_amx(ctxt,tslv)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)
  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Numer of levels of aggr. hierarchy: ",i12)') prec%get_nlevs()
    write(psb_out_unit,'("Time to build aggr. hierarchy     : ",es12.5)')  thier
    write(psb_out_unit,'("Time to build smoothers           : ",es12.5)')  tprec
    write(psb_out_unit,'("Total preconditioner time         : ",es12.5)')  tprec+thier
    write(psb_out_unit,'("Time to solve system              : ",es12.5)')  tslv
    write(psb_out_unit,'("Time per iteration                : ",es12.5)')  tslv/iter
    write(psb_out_unit,'("Number of iterations              : ",i0)')      iter
    write(psb_out_unit,'("Convergence indicator on exit     : ",es12.5)')  err
    write(psb_out_unit,'("Info  on exit                     : ",i0)')      info
    write(psb_out_unit,'("Total memory occupation for      A: ",i12)') amatsize
    write(psb_out_unit,'("Storage format for               A: ",a)')   trim(a%get_fmt())
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)') descsize
    write(psb_out_unit,'("Storage format for          DESC_A: ",a)')   trim(desc_a%get_fmt())
    write(psb_out_unit,'("Total memory occupation for   PREC: ",i12)') precsize
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
  subroutine  get_parms(ctxt,kmethd,afmt,idim,istopc,itmax,itrace,irst,eps)
    
    type(psb_ctxt_type) :: ctxt
    character(len=*)    :: kmethd, afmt
    integer(psb_ipk_)   :: idim, istopc,itmax,itrace,irst
    integer(psb_ipk_)   :: np, iam, info, inp_unit
    real(psb_dpk_)      :: eps
    character(len=20)   :: buffer
    character(len=1024)   :: filename

    call psb_info(ctxt, iam, np)


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
      call read_data(kmethd,inp_unit)
      call read_data(afmt,inp_unit)
      call read_data(idim,inp_unit)
      call read_data(istopc,inp_unit)
      call read_data(itmax,inp_unit)
      call read_data(itrace,inp_unit)
      call read_data(irst,inp_unit)
      call read_data(eps,inp_unit)

      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    ! broadcast parameters to all processors
    call psb_bcast(ctxt,kmethd)
    call psb_bcast(ctxt,afmt)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,istopc)
    call psb_bcast(ctxt,itmax)
    call psb_bcast(ctxt,itrace)
    call psb_bcast(ctxt,irst)
    call psb_bcast(ctxt,eps)

    if (iam == psb_root_) then 
      write(psb_out_unit,'("Solving matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)') np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'("Preconditioner       : ",a)') 'ML-TLU'
      write(psb_out_unit,'("Iterative method     : ",a)') kmethd
      write(psb_out_unit,'(" ")')
    endif

    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer(psb_ipk_) :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  amg_d_pde3d methd prec dim &
         &[istop itmax itrace]'  
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl' 
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '  
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage

end program amg_d_pde3d

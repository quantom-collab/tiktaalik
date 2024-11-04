! matevo.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created June 6, 2024

module matevo
  use pixelation, only: initialize_lagrange_weights
  use alpha_qcd,  only: get_alpha_QCD, get_neff
  use convolution
  use kernels_common
  use kernels_lo

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: pi = acos(-1.0_dp)

  integer, parameter, private :: nfl_min = 3
  integer, parameter, private :: nfl_max = 5

  ! Used for caching
  integer  :: nx_cache = 0, nxi_cache = 0, nQ2_cache = 0

  real(dp), allocatable, dimension(:) :: xi_cache

  ! Zero kernel is just nx by nx
  real(dp), allocatable, dimension(:,:) :: K_zero, K_zero_2

  ! In kernels, indices are for (x,y,xi,nfl)
  ! For singlet/gluon, the first two indices go up to 2*nx,
  ! with the first nx values being singlet quark and the last nx being gluon.
  real(dp), allocatable, dimension(:,:,:,:) :: K_NS_0, KV_SG_0, KA_SG_0

  ! In evolution matrices, indices are for (x,y,xi,Q2)
  ! For singlet/gluon, the first two indices go up to 2*nx,
  ! with the first nx values being singlet quark and the last nx being gluon.
  real(dp), allocatable, dimension(:,:,:,:) :: MV_NS, MV_SG, MA_SG

  public :: make_kernels, make_matrices, &
      & evomat_V_NS, evomat_V_SG, evomat_A_NS, evomat_A_SG, &
      & kernel_V_QQ, kernel_V_QG, kernel_V_GQ, kernel_V_GG, &
      & get_nx, get_nxi, get_nQ2

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to build the kernels and evolution matrices.
    ! IMPORTANT:
    ! - make_kernels **MUST** be called before any matrices can be made!
    !   and also before any routines to grab kernels are called.
    ! - make_matrices **MUST** be called before any routines to grab the
    !   evolution matrices are called.

    subroutine make_kernels(nx, nxi, xi_array)
        ! Public method to initialize the kernel matrices
        ! This **MUST** be called before any evolution matrix routines!
        integer,  intent(in) :: nx, nxi
        real(dp), intent(in) :: xi_array(nxi)
        ! Everything is deallocated first as a safety measure
        call deallocate_all()
        ! Initialize Lagrange weights
        call initialize_lagrange_weights(nx, 6)
        ! Initialize the 2D grids
        call initialize_2D(nx, nxi, xi_array)
        ! Initialize the kernel matrices
        call make_kernels_NS_0(nx, nxi)
        call make_kernels_SG_0(nx, nxi)
    end subroutine make_kernels

    subroutine make_matrices(nQ2, Q2_array)
        ! Public method to initialize the evolution matrices
        ! This **MUST** be called before anything to grab the evolution matrices!
        integer,  intent(in) :: nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        ! Keep track of the number of Q2 points
        nQ2_cache = nQ2
        ! Make evolution matrices
        call make_evomat_NS(nx_cache, nxi_cache, nQ2, Q2_array)
        call make_evomat_V_SG(nx_cache, nxi_cache, nQ2, Q2_array)
        call make_evomat_A_SG(nx_cache, nxi_cache, nQ2, Q2_array)
    end subroutine make_matrices

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Methods to tell Python what the currently cached grid sizes are

    function get_nx() result(nx)
        integer :: nx
        nx = nx_cache
    end function get_nx

    function get_nxi() result(nxi)
        integer :: nxi
        nxi = nxi_cache
    end function get_nxi

    function get_nQ2() result(nQ2)
        integer :: nQ2
        nQ2 = nQ2_cache
    end function get_nQ2

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to return kenrel matrices

    function kernel_V_QQ(nx, nxi, nfl) result(K)
        integer,  intent(in) :: nx, nxi, nfl
        real(dp), dimension(nx,nx,nxi) :: K
        K = K_NS_0(:,:,:,nfl)
    end function kernel_V_QQ

    function kernel_V_QG(nx, nxi, nfl) result(K)
        integer,  intent(in) :: nx, nxi, nfl
        real(dp), dimension(nx,nx,nxi) :: K
        K = KV_SG_0(1:nx,nx+1:2*nx,:,nfl)
    end function kernel_V_QG

    function kernel_V_GQ(nx, nxi, nfl) result(K)
        integer,  intent(in) :: nx, nxi, nfl
        real(dp), dimension(nx,nx,nxi) :: K
        K = KV_SG_0(nx+1:2*nx,1:nx,:,nfl)
    end function kernel_V_GQ

    function kernel_V_GG(nx, nxi, nfl) result(K)
        integer,  intent(in) :: nx, nxi, nfl
        real(dp), dimension(nx,nx,nxi) :: K
        K = KV_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl)
    end function kernel_V_GG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to return evolution matrices

    function evomat_V_NS(nx, nxi, nQ2) result(M)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), dimension(nx, nx, nxi, nQ2) :: M
        M = MV_NS
    end function evomat_V_NS

    function evomat_V_SG(nx, nxi, nQ2) result(M)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), dimension(2*nx, 2*nx, nxi, nQ2) :: M
        M = MV_SG
    end function evomat_V_SG

    function evomat_A_NS(nx, nxi, nQ2) result(M)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), dimension(nx, nx, nxi, nQ2) :: M
        M = MV_NS
    end function evomat_A_NS

    function evomat_A_SG(nx, nxi, nQ2) result(M)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), dimension(2*nx, 2*nx, nxi, nQ2) :: M
        M = MA_SG
    end function evomat_A_SG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Methods to make kernel matrices

    subroutine make_kernels_NS_0(nx, nxi)
        integer, intent(in) :: nx, nxi
        integer :: ix, iy, iz!, nfl
        if(allocated(K_NS_0)) deallocate(K_NS_0)
        allocate(K_NS_0(nx,nx,nxi,nfl_min:nfl_max))
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              ! At leading order, no nfl dependence in QQ kernel.
              K_NS_0(ix,iy,iz,:) = pixel_conv(zero_func, K0_qq_pls, K0_qq_cst, xi_cache(iz), nx, ix, iy)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
    end subroutine make_kernels_NS_0

    subroutine make_kernels_SG_0(nx, nxi)
        integer, intent(in) :: nx, nxi
        integer :: ix, iy, iz, nfl
        real(dp), dimension(:,:,:), allocatable :: qq_nfl_0, qG_nfl_1, Gq_nfl_0, GG_nfl_0, GG_nfl_1
        if(allocated(KA_SG_0)) deallocate(KA_SG_0)
        if(allocated(KV_SG_0)) deallocate(KV_SG_0)
        allocate(KA_SG_0(2*nx,2*nx,nxi,nfl_min:nfl_max))
        allocate(KV_SG_0(2*nx,2*nx,nxi,nfl_min:nfl_max))
        ! Temporary sub-matrix arrays
        allocate(qq_nfl_0(nx,nx,nxi))
        allocate(qG_nfl_1(nx,nx,nxi))
        allocate(Gq_nfl_0(nx,nx,nxi))
        allocate(GG_nfl_0(nx,nx,nxi))
        allocate(GG_nfl_1(nx,nx,nxi))
        KA_SG_0 = 0.0_dp
        KV_SG_0 = 0.0_dp
        ! First, build up the A-type kernel
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              qG_nfl_1(ix,iy,iz) = pixel_conv(KA0_qG_reg, zero_func,  zero_func,  xi_cache(iz), nx, ix, iy)
              Gq_nfl_0(ix,iy,iz) = pixel_conv(KA0_Gq_reg, zero_func,  zero_func,  xi_cache(iz), nx, ix, iy)
              GG_nfl_0(ix,iy,iz) = pixel_conv(zero_func,  KA0_GG_pls, KA0_GG_cst, xi_cache(iz), nx, ix, iy)
              GG_nfl_1(ix,iy,iz) = pixel_conv(zero_func,  zero_func,  KA0_GG_nfl, xi_cache(iz), nx, ix, iy)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        ! For the QQ block, just copy over the NS kernel, since it's the same at LO.
        ! The NS kernel matrix is built before this method is called so this is safe.
        KA_SG_0(1:nx,1:nx,:,:) = K_NS_0(:,:,:,:)
        do nfl=nfl_min, nfl_max, 1
          !!KA_SG_0(1:nx,     1:nx,     :,nfl) = qq_nfl_0(:,:,:)
          KA_SG_0(1:nx,     nx+1:2*nx,:,nfl) = real(nfl)*qG_nfl_1(:,:,:)
          KA_SG_0(nx+1:2*nx,1:nx,     :,nfl) = Gq_nfl_0(:,:,:)
          KA_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl) = GG_nfl_0(:,:,:) + real(nfl)*GG_nfl_1(:,:,:)
        end do
        ! Next, the extra pieces for the V-type kernel
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              qG_nfl_1(ix,iy,iz) = pixel_conv(KVmA0_qG_reg, zero_func, zero_func, xi_cache(iz), nx, ix, iy)
              Gq_nfl_0(ix,iy,iz) = pixel_conv(KVmA0_Gq_reg, zero_func, zero_func, xi_cache(iz), nx, ix, iy)
              GG_nfl_0(ix,iy,iz) = pixel_conv(KVmA0_GG_reg, zero_func, zero_func, xi_cache(iz), nx, ix, iy)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        KV_SG_0(1:nx,1:nx,:,:) = 0.0_dp
        do nfl=nfl_min, nfl_max, 1
          KV_SG_0(1:nx,     nx+1:2*nx,:,nfl) = real(nfl)*qG_nfl_1(:,:,:)
          KV_SG_0(nx+1:2*nx,1:nx,     :,nfl) = Gq_nfl_0(:,:,:)
          KV_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl) = GG_nfl_0(:,:,:)
        end do
        ! And add the extra pieces to the A-type kernel
        KV_SG_0 = KV_SG_0 + KA_SG_0
        ! And we're done
        deallocate(qq_nfl_0)
        deallocate(qG_nfl_1)
        deallocate(Gq_nfl_0)
        deallocate(GG_nfl_0)
        deallocate(GG_nfl_1)
    end subroutine make_kernels_SG_0

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Methods to make evolution matrices

    subroutine make_evomat_NS(nx, nxi, nQ2, Q2_array)
        ! Just one non-singlet evolution matrix at leading order.
        ! I'll need to break this into multiple routines,
        ! for both V-type and A-type, but also plus-type and minus-type,
        ! after going to NLO. But I'll burn that bridge when I come to it.
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        real(dp), dimension(nx,nx) :: idnx
        integer :: ix, ixi, iQ2, nfl
        if(allocated(MV_NS)) deallocate(MV_NS)
        allocate(MV_NS(nx,nx,nxi,nQ2))
        ! Build an identity matrix
        idnx = 0.0_dp
        do ix=1, nx, 1
          idnx(ix,ix) = 1.0_dp
        end do
        ! Identity for no evolution at initial scale
        MV_NS = 0.0_dp
        do ixi=1, nxi, 1
          MV_NS(:,:,ixi,1) = idnx
        end do
        ! Build evolution matrices for all other scales
        do iQ2=2, nQ2, 1
          !$OMP PARALLEL DO
          do ixi=1, nxi, 1
            ! First, an evolution matrix for prior Q2 step to current Q2 step
            nfl = get_neff(Q2_array(iQ2-1))
            MV_NS(:,:,ixi,iQ2) = idnx + &
                & rk4_NS(nx, nxi, Q2_array(iQ2-1), Q2_array(iQ2), K_NS_0(:,:,ixi,nfl), K_zero(:,:))
            ! Then matrix multiplication to turn into matrix from initial Q2 to current Q2
            MV_NS(:,:,ixi,iQ2) = matmul(MV_NS(:,:,ixi,iQ2), MV_NS(:,:,ixi,iQ2-1))
          end do
          !$OMP END PARALLEL DO
        end do
    end subroutine make_evomat_NS

    subroutine make_evomat_V_SG(nx, nxi, nQ2, Q2_array)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        real(dp), dimension(2*nx,2*nx) :: id2nx
        integer :: ix, ixi, iQ2, nfl
        if(allocated(MV_SG)) deallocate(MV_SG)
        allocate(MV_SG(2*nx,2*nx,nxi,nQ2))
        ! Build an identity matrix
        id2nx = 0.0_dp
        do ix=1, 2*nx, 1
          id2nx(ix,ix) = 1.0_dp
        end do
        ! Identity for no evolution at initial scale
        MV_SG = 0.0_dp
        do ixi=1, nxi, 1
          MV_SG(:,:,ixi,1) = id2nx
        end do
        ! Build evolution matrices for all other scales
        do iQ2=2, nQ2, 1
          !$OMP PARALLEL DO
          do ixi=1, nxi, 1
            ! First, an evolution matrix for prior Q2 step to current Q2 step
            nfl = get_neff(Q2_array(iQ2-1))
            MV_SG(:,:,ixi,iQ2) = id2nx + &
                & rk4_SG(nx, nxi, Q2_array(iQ2-1), Q2_array(iQ2), &
                & KV_SG_0(:,:,ixi,nfl), K_zero_2(:,:))
            ! Then matrix multiplication to turn into matrix from initial Q2 to current Q2
            MV_SG(:,:,ixi,iQ2) = matmul(MV_SG(:,:,ixi,iQ2), MV_SG(:,:,ixi,iQ2-1))
          end do
          !$OMP END PARALLEL DO
        end do
    end subroutine make_evomat_V_SG

    subroutine make_evomat_A_SG(nx, nxi, nQ2, Q2_array)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        real(dp), dimension(2*nx,2*nx) :: id2nx
        integer :: ix, ixi, iQ2, nfl
        if(allocated(MA_SG)) deallocate(MA_SG)
        allocate(MA_SG(2*nx,2*nx,nxi,nQ2))
        ! Build an identity matrix
        id2nx = 0.0_dp
        do ix=1, 2*nx, 1
          id2nx(ix,ix) = 1.0_dp
        end do
        ! Identity for no evolution at initial scale
        MA_SG = 0.0_dp
        do ixi=1, nxi, 1
          MA_SG(:,:,ixi,1) = id2nx
        end do
        ! Build evolution matrices for all other scales
        do iQ2=2, nQ2, 1
          !$OMP PARALLEL DO
          do ixi=1, nxi, 1
            ! First, an evolution matrix for prior Q2 step to current Q2 step
            nfl = get_neff(Q2_array(iQ2-1))
            MA_SG(:,:,ixi,iQ2) = id2nx + &
                & rk4_SG(nx, nxi, Q2_array(iQ2-1), Q2_array(iQ2), &
                & KA_SG_0(:,:,ixi,nfl), K_zero_2(:,:))
            ! Then matrix multiplication to turn into matrix from initial Q2 to current Q2
            MA_SG(:,:,ixi,iQ2) = matmul(MA_SG(:,:,ixi,iQ2), MA_SG(:,:,ixi,iQ2-1))
          end do
          !$OMP END PARALLEL DO
        end do
    end subroutine make_evomat_A_SG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! RK4

    function rk4_NS(nx, nxi, Q2i, Q2f, K0, K1) result(K)
    !subroutine rk4_NS(nx, nxi, Q2i, Q2f, K0, K1, K)
        integer,  intent(in) :: nx, nxi
        real(dp), intent(in) :: Q2i, Q2f, K0(nx,nx), K1(nx,nx)
        !real(dp), intent(out) :: K(nx,nx)
        real(dp) :: K(nx,nx)
        !
        !real(dp), dimension(nx,nx) :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp), dimension(:,:), allocatable :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp) :: h, Q2m, a_ini, a_mid, a_end
        allocate(W1(nx,nx), W2(nx,nx), W3(nx,nx), W4(nx,nx), M_ini(nx,nx), M_mid(nx,nx), M_end(nx,nx))
        h = log(Q2f/Q2i)
        Q2m = sqrt(Q2i*Q2f) ! geometric mean for midpoint
        !
        a_ini = get_alpha_QCD(Q2i) / (2.*pi)
        a_mid = get_alpha_QCD(Q2m) / (2.*pi)
        a_end = get_alpha_QCD(Q2f) / (2.*pi)
        ! The three matrices over each subinterval
        M_ini = a_ini*K0 + a_ini**2*K1
        M_mid = a_mid*K0 + a_mid**2*K1
        M_end = a_end*K0 + a_end**2*K1
        ! The RK4 k-values
        W1 = M_ini
        W2 = M_mid + 0.5*h*matmul(M_mid,W1)
        W3 = M_mid + 0.5*h*matmul(M_mid,W2)
        W4 = M_end + h*matmul(M_end,W3)
        ! And the final formula for the shift
        K = h/6.*(W1 + 2.*W2 + 2.*W3 + W4)
        deallocate(W1, W2, W3, W4, M_ini, M_mid, M_end)
    end function rk4_NS
    !end subroutine rk4_NS

    function rk4_SG(nx, nxi, Q2i, Q2f, K0, K1) result(K)
    !subroutine rk4_SG(nx, nxi, Q2i, Q2f, K0, K1, K)
        integer,  intent(in) :: nx, nxi
        real(dp), intent(in) :: Q2i, Q2f, K0(2*nx,2*nx), K1(2*nx,2*nx)
        !real(dp), intent(out) :: :: K(2*nx,2*nx)
        real(dp) :: K(2*nx,2*nx)
        !
        !real(dp), dimension(2*nx,2*nx) :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp), dimension(:,:), allocatable :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp) :: h, Q2m, a_ini, a_mid, a_end
        allocate(W1(2*nx,2*nx), W2(2*nx,2*nx), W3(2*nx,2*nx), W4(2*nx,2*nx), M_ini(2*nx,2*nx), M_mid(2*nx,2*nx), M_end(2*nx,2*nx))
        h = log(Q2f/Q2i)
        Q2m = sqrt(Q2i*Q2f) ! geometric mean for midpoint
        !
        a_ini = get_alpha_QCD(Q2i) / (2.*pi)
        a_mid = get_alpha_QCD(Q2m) / (2.*pi)
        a_end = get_alpha_QCD(Q2f) / (2.*pi)
        ! The three matrices over each subinterval
        M_ini = a_ini*K0 + a_ini**2*K1
        M_mid = a_mid*K0 + a_mid**2*K1
        M_end = a_end*K0 + a_end**2*K1
        ! The RK4 k-values
        W1 = M_ini
        W2 = M_mid + 0.5*h*matmul(M_mid,W1)
        W3 = M_mid + 0.5*h*matmul(M_mid,W2)
        W4 = M_end + h*matmul(M_end,W3)
        ! And the final formula for the shift
        K = h/6.*(W1 + 2.*W2 + 2.*W3 + W4)
        deallocate(W1, W2, W3, W4, M_ini, M_mid, M_end)
    end function rk4_SG
    !end subroutine rk4_SG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Methods to initialize some cached arrays

    subroutine initialize_2D(nx, nxi, xi_grid)
        integer, intent(in)  :: nx, nxi
        real(dp), intent(in) :: xi_grid(nxi)
        nx_cache  = nx
        nxi_cache = nxi
        if(allocated(xi_cache)) deallocate(xi_cache)
        allocate(xi_cache(nxi))
        xi_cache = xi_grid
        ! To help with singular behavior ...
        if(xi_cache(1) < 1e-6_dp) xi_cache(1) = 1e-6_dp
        if(xi_cache(nxi) > 0.999999_dp) xi_cache(nxi) = 0.999999_dp
        if(allocated(K_zero)) deallocate(K_zero)
        allocate(K_zero(nx,nx))
        K_zero = 0.0_dp
        if(allocated(K_zero_2)) deallocate(K_zero_2)
        allocate(K_zero_2(2*nx,2*nx))
        K_zero_2 = 0.0_dp
    end subroutine initialize_2D

    subroutine deallocate_all()
        if(allocated(K_NS_0))  deallocate(K_NS_0)
        if(allocated(KV_SG_0)) deallocate(KV_SG_0)
        if(allocated(MV_NS))   deallocate(MV_NS)
        if(allocated(MV_SG))   deallocate(MV_SG)
        if(allocated(KA_SG_0)) deallocate(KA_SG_0)
        if(allocated(MA_SG))   deallocate(MA_SG)
        nQ2_cache = 0
    end subroutine deallocate_all

end module matevo

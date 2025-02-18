! matrices.f90
!
! by Adma Freese
! part of the package tiktaalik
!
! wrappers for f2py to access

module dummy
  use gridspace
  use matevo
  use pixelation
  use specfun ! TODO REMOVE BEFORE RELEASE

  implicit none
  public

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Pixelspace

    subroutine pixelspace_wrap(nx, xi, grid_type, x)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, grid_type
        real(dp), intent(in)  :: xi
        real(dp), intent(out) :: x(nx)
        integer :: ix
        do ix=1, nx, 1
          x(ix) = push_forward(real(2*ix-1)/real(nx) - 1.0_dp, xi, nx, grid_type)
        end do
    end subroutine pixelspace_wrap

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialization routines, **MUST** be called first!

    subroutine make_kernels_wrap(nx, nxi, xi_array, grid_type)
        ! Initializes kernel matrices, using a particular xi array.
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nx, nxi, grid_type
        real(dp), intent(in) :: xi_array(nxi)
        call make_kernels(nx, nxi, xi_array, grid_type)
    end subroutine make_kernels_wrap

    subroutine make_matrices_wrap(nQ2, Q2_array, l_nlo)
        ! Initializes evolution matrices, using a particular Q2 array.
        ! The kernels must have already been initialized.
        ! The nx and nxi used here must be consistent with the nx and nxi
        ! that the kernels were initialized with.
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        logical,  intent(in) :: l_nlo
        call make_matrices(nQ2, Q2_array, l_nlo)
    end subroutine make_matrices_wrap

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Routines to pass chached grid sizes to Python (to ease user burden)

    subroutine get_nx_wrap(nx)
        integer, intent(out) :: nx
        nx = get_nx()
    end subroutine get_nx_wrap

    subroutine get_nxi_wrap(nxi)
        integer, intent(out) :: nxi
        nxi = get_nxi()
    end subroutine get_nxi_wrap

    subroutine get_nQ2_wrap(nQ2)
        integer, intent(out) :: nQ2
        nQ2 = get_nQ2()
    end subroutine get_nQ2_wrap

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Evolution matrices

    subroutine evomatrix_vns_wrap(nx, nxi, nQ2, nstype, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2, nstype
        real(dp), intent(out) :: M(nx,nx,nxi,nQ2)
        !
        M = evomat_V_NS(nx, nxi, nQ2, nstype)
    end subroutine evomatrix_vns_wrap

    subroutine evomatrix_vsg_wrap(nx, nxi, nQ2, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        real(dp), intent(out) :: M(2*nx,2*nx,nxi,nQ2)
        !
        M = evomat_V_SG(nx, nxi, nQ2)
    end subroutine evomatrix_vsg_wrap

    subroutine evomatrix_ans_wrap(nx, nxi, nQ2, nstype, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2, nstype
        real(dp), intent(out) :: M(nx,nx,nxi,nQ2)
        !
        M = evomat_A_NS(nx, nxi, nQ2, nstype)
    end subroutine evomatrix_ans_wrap

    subroutine evomatrix_asg_wrap(nx, nxi, nQ2, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        real(dp), intent(out) :: M(2*nx,2*nx,nxi,nQ2)
        !
        M = evomat_A_SG(nx, nxi, nQ2)
    end subroutine evomatrix_asg_wrap

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Kernel matrices

    subroutine evokernel_vqq_wrap(Q2, nx, nxi, nfl, l_nlo, i_ns_type, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl, i_ns_type
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_qq(Q2, nx, nxi, nfl, l_nlo, i_ns_type)
    end subroutine evokernel_vqq_wrap

    subroutine evokernel_vqg_wrap(Q2, nx, nxi, nfl, l_nlo, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_qg(Q2, nx, nxi, nfl, l_nlo)
    end subroutine evokernel_vqg_wrap

    subroutine evokernel_vgq_wrap(Q2, nx, nxi, nfl, l_nlo, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_gq(Q2, nx, nxi, nfl, l_nlo)
    end subroutine evokernel_vgq_wrap

    subroutine evokernel_vgg_wrap(Q2, nx, nxi, nfl, l_nlo, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_gg(Q2, nx, nxi, nfl, l_nlo)
    end subroutine evokernel_vgg_wrap

    subroutine evokernel_aqq_wrap(Q2, nx, nxi, nfl, l_nlo, i_ns_type, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl, i_ns_type
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_A_qq(Q2, nx, nxi, nfl, l_nlo, i_ns_type)
    end subroutine evokernel_aqq_wrap

    subroutine evokernel_aqg_wrap(Q2, nx, nxi, nfl, l_nlo, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_A_qg(Q2, nx, nxi, nfl, l_nlo)
    end subroutine evokernel_aqg_wrap

    subroutine evokernel_agq_wrap(Q2, nx, nxi, nfl, l_nlo, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_A_gq(Q2, nx, nxi, nfl, l_nlo)
    end subroutine evokernel_agq_wrap

    subroutine evokernel_agg_wrap(Q2, nx, nxi, nfl, l_nlo, K)
        integer,  parameter   :: dp = kind(1d0)
        real(dp), intent(in)  :: Q2
        integer,  intent(in)  :: nx, nxi, nfl
        logical,  intent(in)  :: l_nlo
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_A_gg(Q2, nx, nxi, nfl, l_nlo)
    end subroutine evokernel_agg_wrap

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! TESTING AREA (TODO REMOVE BEFORE RELEASE)

    subroutine dilog_wrap(nx, x, y)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx
        real(dp), intent(in)  :: x(nx)
        real(dp), intent(out) :: y(nx)
        !
        integer :: i
        do i=1, nx, 1
          y(i) = dilog(x(i))
        end do
    end subroutine dilog_wrap

end module dummy

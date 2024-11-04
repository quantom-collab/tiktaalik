! matrices.f90
!
! by Adma Freese
! part of the package tiktaalik
!
! wrappers for f2py to access

module dummy
  use pixelation
  use matevo

  implicit none
  public

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialization routines, **MUST** be called first!

    subroutine make_kernels_wrap(nx, nxi, xi_array)
        ! Initializes kernel matrices, using a particular xi array.
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nx, nxi
        real(dp), intent(in) :: xi_array(nxi)
        call make_kernels(nx, nxi, xi_array)
    end subroutine make_kernels_wrap

    subroutine make_matrices_wrap(nQ2, Q2_array)
        ! Initializes evolution matrices, using a particular Q2 array.
        ! The kernels must have already been initialized.
        ! The nx and nxi used here must be consistent with the nx and nxi
        ! that the kernels were initialized with.
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        call make_matrices(nQ2, Q2_array)
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

    subroutine evomatrix_vns_wrap(nx, nxi, nQ2, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        real(dp), intent(out) :: M(nx,nx,nxi,nQ2)
        !
        M = evomat_V_NS(nx, nxi, nQ2)
    end subroutine evomatrix_vns_wrap

    subroutine evomatrix_vsg_wrap(nx, nxi, nQ2, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        real(dp), intent(out) :: M(2*nx,2*nx,nxi,nQ2)
        !
        M = evomat_V_SG(nx, nxi, nQ2)
    end subroutine evomatrix_vsg_wrap

    subroutine evomatrix_ans_wrap(nx, nxi, nQ2, M)
        ! QQ, helicity-independent
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        real(dp), intent(out) :: M(nx,nx,nxi,nQ2)
        !
        M = evomat_A_NS(nx, nxi, nQ2)
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

    subroutine evokernel_vqq_wrap(nx, nxi, nfl, K)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nfl
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_qq(nx, nxi, nfl)
    end subroutine evokernel_vqq_wrap

    subroutine evokernel_vqg_wrap(nx, nxi, nfl, K)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nfl
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_qg(nx, nxi, nfl)
    end subroutine evokernel_vqg_wrap

    subroutine evokernel_vgq_wrap(nx, nxi, nfl, K)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nfl
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_gq(nx, nxi, nfl)
    end subroutine evokernel_vgq_wrap

    subroutine evokernel_vgg_wrap(nx, nxi, nfl, K)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nfl
        real(dp), intent(out) :: K(nx,nx,nxi)
        !
        K = kernel_V_gg(nx, nxi, nfl)
    end subroutine evokernel_vgg_wrap

end module dummy

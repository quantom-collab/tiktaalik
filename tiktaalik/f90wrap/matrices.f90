! matrices.f90
!
! by Adma Freese
! part of the package tiktaalik
!
! wrappers for f2py to access

module dummy
  use gridspace
  use matrices_common
  use matrices_evolution
  use matrices_mellin
  use matrices_wilson
  use pixelation

  implicit none
  public

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Interpixel

    subroutine interpixel_wrap(n_pixels, i_pixel, nx, x, xi, grid_type, y)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: n_pixels, i_pixel, nx, grid_type
        real(dp), intent(in)  :: x(nx), xi
        real(dp), intent(out) :: y(nx)
        !
        integer :: ix
        do ix=1, nx, 1
          y(ix) = interpixel(n_pixels, i_pixel, x(ix), xi, grid_type)
        end do
    end subroutine interpixel_wrap

    ! TODO : deprecate
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
    ! Initialization routines

    subroutine initialize_x_xi_wrap(nx, nxi, xi_grid, grid_type, lagrange_order)
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nx, nxi, grid_type, lagrange_order
        real(dp), intent(in) :: xi_grid(nxi)
        call initialize_x_xi(nx, nxi, xi_grid, grid_type, lagrange_order)
    end subroutine initialize_x_xi_wrap

    subroutine initialize_Q2_wrap(nQ2, Q2_array)
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        call initialize_Q2(nQ2, Q2_array)
    end subroutine initialize_Q2_wrap

    subroutine make_kernels_wrap()
        ! Initializes kernel matrices
        call make_kernels()
    end subroutine make_kernels_wrap

    subroutine make_matrices_wrap(nQ2, Q2_array, l_nlo)
        ! Initializes evolution matrices, using a particular Q2 array.
        ! The kernels must have already been initialized.
        integer,  parameter  :: dp = kind(1d0)
        integer,  intent(in) :: nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        logical,  intent(in) :: l_nlo
        call make_evolution_matrices(nQ2, l_nlo)
    end subroutine make_matrices_wrap

    subroutine make_wilson_wrap(nQ2)
        ! Initializes Wilson coefficient matrices.
        integer,  intent(in) :: nQ2
        call make_wilson_matrices(nQ2)
    end subroutine make_wilson_wrap

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

    subroutine get_x_wrap(nx, nxi, xx)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(out) :: xx(nx,nxi)
        xx = get_x(nx, nxi)
    end subroutine get_x_wrap

    subroutine get_xi_wrap(nxi, xi)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nxi
        real(dp), intent(out) :: xi(nxi)
        xi = get_xi(nxi)
    end subroutine get_xi_wrap

    subroutine get_Q2_wrap(nQ2, Q2)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nQ2
        real(dp), intent(out) :: Q2(nQ2)
        Q2 = get_Q2(nQ2)
    end subroutine get_Q2_wrap

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
    ! Wilson coefficient matrices

    subroutine dvcs_cq_wrap(nx, nxi, nQ2, l_nlo, C)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        logical,  intent(in)  :: l_nlo
        complex(dp), intent(out) :: C(nxi, nx, nQ2)
        !
        C = Cq_dvcs(nxi, nx, nQ2, l_nlo)
    end subroutine dvcs_cq_wrap

    subroutine dvcs_cg_wrap(nx, nxi, nQ2, l_nlo, C)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        logical,  intent(in)  :: l_nlo
        complex(dp), intent(out) :: C(nxi, nx, nQ2)
        !
        C = CG_dvcs(nxi, nx, nQ2, l_nlo)
    end subroutine dvcs_cg_wrap

    subroutine dvcs_ctilq_wrap(nx, nxi, nQ2, l_nlo, C)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        logical,  intent(in)  :: l_nlo
        complex(dp), intent(out) :: C(nxi, nx, nQ2)
        !
        C = Ctilq_dvcs(nxi, nx, nQ2, l_nlo)
    end subroutine dvcs_ctilq_wrap

    subroutine dvcs_ctilg_wrap(nx, nxi, nQ2, l_nlo, C)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nQ2
        logical,  intent(in)  :: l_nlo
        complex(dp), intent(out) :: C(nxi, nx, nQ2)
        !
        C = CtilG_dvcs(nxi, nx, nQ2, l_nlo)
    end subroutine dvcs_ctilg_wrap
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Mellin matrices

    subroutine mellin_wrap(nx, nxi, s, M)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, s
        real(dp), intent(out) :: M(nxi, nx)
        !
        call make_mellin_matrix(nx, nxi, s, M)
    end subroutine mellin_wrap

end module dummy

! matrices_mellin.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created May 21, 2026
! To get Mellin moments for obtaining generalized form factors,
! and for testing polynomality

module matrices_mellin
  use matrices_common
  use kernels_common,  only: zero_func
  use wilson_dvcs,     only: pixel_coef

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  public :: make_mellin_matrix

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods

    subroutine make_mellin_matrix(nx, nxi, s, M)
        integer,  intent(in)  :: nx, nxi, s
        real(dp), intent(out) :: M(nxi,nx)
        real(dp) :: xi(nxi)
        integer :: grid_type, ix, iz
        xi = get_xi(nxi)
        grid_type = get_grid_type()
        do iz=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          M(iz,ix) = pixel_coef(mellin_kernel, zero_func, zero_func, xi(iz), nx, ix, grid_type)
        end do
        !$OMP END PARALLEL DO
        end do
        return
        contains
          function mellin_kernel(x) result(y)
              real(dp), intent(in) :: x
              real(dp) :: y
              y = x**(s-1)
          end function mellin_kernel
    end subroutine make_mellin_matrix

end module matrices_mellin

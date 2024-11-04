! convolution.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! this comment added April 4, 2024, but parts of the the code were written over
! the course of 2024.

module convolution
  use integration, only: integrate
  use pixelation,  only: interpixel, linspace

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  public :: pixel_conv

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Pixelated shifts

    function pixel_conv(Kreg, Kpls, Kcst, xi, N, i, j) result(shift)
        real(dp), external   :: Kreg, Kpls, Kcst
        real(dp), intent(in) :: xi
        integer,  intent(in) :: N, i, j
        real(dp) :: shift
        !
        shift = 0.0_dp
        shift = shift + pixel_conv_reg(Kreg, xi, N, i, j)
        shift = shift + pixel_conv_pls(Kpls, xi, N, i, j)
        shift = shift + pixel_conv_cst(Kcst, xi, N, i, j)
    end function pixel_conv

    function pixel_conv_reg(kernel, xi, n_pixels, i, j) result(shift)
        real(dp), external   :: kernel
        real(dp), intent(in) :: xi
        integer,  intent(in) :: n_pixels, i, j
        real(dp) :: shift
        !
        real(dp) :: x
        x = real(2*i-1)/real(n_pixels) - 1.
        shift = integrate(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              real(dp) :: fy
              fy = interpixel(n_pixels, j, y)
              intd = kernel(x,y,xi)*fy
          end function integrand
    end function pixel_conv_reg

    function pixel_conv_pls(kernel, xi, n_pixels, i, j) result(shift)
        real(dp), external   :: kernel
        real(dp), intent(in) :: xi
        integer,  intent(in) :: n_pixels, i, j
        real(dp) :: shift
        !
        real(dp) :: x
        x = real(2*i-1)/real(n_pixels) - 1.
        shift = integrate(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              real(dp) :: fy, fx
              fy = interpixel(n_pixels, j, y)
              fx = interpixel(n_pixels, j, x)
              intd = kernel(x,y,xi)*(fy-fx)
          end function integrand
    end function pixel_conv_pls

    function pixel_conv_cst(kernel, xi, n_pixels, i, j) result(shift)
        real(dp), external   :: kernel
        real(dp), intent(in) :: xi
        integer,  intent(in) :: n_pixels, i, j
        real(dp) :: shift
        !
        real(dp) :: x, fx
        shift = 0.0_dp
        if(i==j) then
          x = real(2*i-1)/real(n_pixels) - 1.
          fx = interpixel(n_pixels, j, x)
          shift = kernel(x,xi) * fx
        endif
    end function pixel_conv_cst

end module convolution

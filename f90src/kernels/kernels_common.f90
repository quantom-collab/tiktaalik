! kernels_common.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! this comment added April 4, 2024, but parts of the the code were written over
! the course of 2024.
!
! All kernels are taken from:
!   Belitsky, Fruend and Mueller
!   Nuclear Physics B 574 (2000) 347-406
!   Belitsky:1999hf
!   arxiv:hep-ph/9912379

module kernels_common
  use constants, only: pi
  use specfun,   only: dilog

  implicit none
  public

  integer,  parameter, private :: dp = kind(1d0)

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Function to serve as zero

    function zero_func() result(z)
        real(dp) :: z
        z = 0.0_dp
    end function zero_func

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Step functions

    function theta_step(X) result(theta)
        real(dp), intent(in) :: X
        real(dp) :: theta
        theta = 0.0_dp
        if(X > 0.0_dp) then
          theta = 1.0_dp
        elseif(X==0.0_dp) then
          theta = 0.5_dp
        endif
    end function theta_step

    function rho_step(X, Y) result(rho)
        real(dp), intent(in) :: X, Y
        real(dp) :: rho
        !
        rho = theta_step(1.-X/Y)*theta_step(X/Y)*sign(1.0_dp,Y)
    end function rho_step

    function erbl_step(X) result(theta)
        ! Check that we're in the ERBL region
        real(dp), intent(in) :: X
        real(dp) :: theta
        !
        theta = theta_step(X)*theta_step(1.0_dp-X)
    end function erbl_step

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Convenience functions for dealing with logs of negative values

    function abslog(x) result(y)
        ! log of absolute value
        ! Convenient for removing imaginary parts
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        y = log(abs(x))
    end function abslog

    function log2(x) result(y)
        ! Square of log
        ! Removes imaginary part, accounts for -pi**2 in real part
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        y = log(abs(x))**2
        if(x < 0.0_dp) y = y - pi**2
    end function log2

    function logprod(x1, x2) result(y)
        ! Product of two logs
        ! Removes imaginary part, accounts for -pi**2 in real part
        real(dp), intent(in) :: x1, x2
        real(dp) :: y
        !
        y = log(abs(x1)) * log(abs(x2))
        if(x1 < 0.0_dp .and. x2 < 0.0_dp) y = y - pi**2
    end function logprod

    function li2(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        y = dilog(x)
    end function li2

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Small pieces of the BFM kernels

    ! QQ

    function QQ_f_a(X,Y) result(f)
        ! Eq. (76)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        f = X/Y
    end function QQ_f_a

    function QQ_f_b(X,Y) result(f)
        ! Eq. (76)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        f = X/Y / (Y-X)
    end function QQ_f_b

    function QQ_f_c(X,Y) result(f)
        ! Eq. (77)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        real(dp) :: Xbar
        Xbar = 1.0_dp - X
        f = X/Y * (2.0_dp*Xbar*Y*(4.0_dp/3.0_dp - abslog(Xbar*Y)) + Y - X)
    end function QQ_f_c

    ! QG

    function QG_f_a(X,Y) result(f)
        ! Eq. (76)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        f = X/Y**2
    end function QG_f_a

    function QG_f_c(X,Y) result(f)
        ! Eq. (78)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        real(dp) :: Xbar
        Xbar = 1.0_dp - X
        f = X/Y**2 * (2.0_dp*Xbar*Y - X)
    end function QG_f_c

    ! GQ

    function GQ_f_a(X,Y) result(f)
        ! Eq. (76)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        f = X**2/Y
    end function GQ_f_a

    function GQ_f_c(X,Y) result(f)
        ! Eq. (78)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        real(dp) :: Xbar, Ybar
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        f = X**2/Y * (2.0_dp*Xbar*Y - Ybar)
    end function GQ_f_c

    function GQ_h_V(X,Y) result(h)
        ! Eq. (174)
        real(dp), intent(in) :: X, Y
        real(dp) :: h
        !
        real(dp) :: Xbar, Ybar, GQfV, GQfVbar
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        ! Eqs. (74), (75)
        GQfV    = GQ_f_a(X   ,Y)    + 2.*GQ_f_c(X   ,Y   )
        GQfVbar = GQ_f_a(Xbar,Ybar) + 2.*GQ_f_c(Xbar,Ybar)
        h = -2.*GQfVbar*logprod(Xbar,Y) - 2.*GQfV*(li2(X) + li2(Ybar))
    end function GQ_h_V

    function GQ_hbar_V(X,Y) result(h)
        ! Eq. (174)
        real(dp), intent(in) :: X, Y
        real(dp) :: h
        !
        real(dp) :: Xbar, Ybar, GQfV, GQfVbar
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        ! Eqs. (74), (75)
        GQfV    = GQ_f_a(X   ,Y)    + 2.*GQ_f_c(X   ,Y   )
        GQfVbar = GQ_f_a(Xbar,Ybar) + 2.*GQ_f_c(Xbar,Ybar)
        ! What's literally written
        h = (GQfVbar + GQfV)*(2.*li2(1.-X/Y) + log2(Y)) &
            & + 2.*GQfV*(li2(Ybar) - logprod(X,Y)) &
            & - 2.*GQfVbar*li2(Xbar)
        ! My guess at what they meant?
        !h = (GQfVbar + GQfV)*(2.*li2(1.-Xbar/Y) + log2(Y)) &
        !    & + 2.*GQfV*(li2(Ybar) - logprod(Xbar,Y)) &
        !    & - 2.*GQfVbar*li2(X)
    end function GQ_hbar_V

    ! GG

    function GG_f_a(X,Y) result(f)
        ! Eq. (76)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        f = X**2/Y**2
    end function GG_f_a

    function GG_f_b(X,Y) result(f)
        ! Eq. (76)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        f = X**2/Y**2 / (Y-X)
    end function GG_f_b

    function GG_f_c(X,Y) result(f)
        ! Eq. (77)
        real(dp), intent(in) :: X, Y
        real(dp) :: f
        !
        real(dp) :: Xbar
        Xbar = 1.0_dp - X
        f = X**2/Y**2 * (2.0_dp*Xbar*Y + Y - X)
    end function GG_f_c

end module kernels_common

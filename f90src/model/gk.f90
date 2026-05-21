! gk.f90
!
! part of the package tiktaalik for GPD evolution
!
! Code added July 30, 2024 by Adam Freese.
!
! This implements a model due to Goloshkokov and Kroll.
! The specific reference I consulted is was:
!   P. Kroll, H. Moutarde, F. Sabatie
!   European Physical Journal C (2013) 73:2278
!   arxiv:1210.6975
!   Kroll:2012sm
! However, the original references for the GK model are:
!   S.V. Goloskokov and P. Kroll,
!   European Physical Journal C (2005) 42:281-301
!   arxiv:hep-ph/0501242
!   Goloskokov:2005sd
! and
!   S.V. Goloskokov and P. Kroll,
!   European Physical Journal C (2008) 53:367-384
!   arxiv:0708.3569
!   Goloskokov:2007nt

module gk
  use gk_analytic, only: H_n1_analytic, H_n2_analytic
  use gk_smol,     only: H_n1_smol,     H_n2_smol
  use gk_zero,     only: H_n1_zero,     H_n2_zero

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: pi = acos(-1.0_dp)
  real(dp), parameter, private :: eps = 1e-12_dp

  ! H-type GPDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! from Table 1 of Kroll:2012sm
  ! using Q = Q0, since the code does evolution

  ! Regge
  real(dp), parameter, private :: a0_val = 0.48_dp
  real(dp), parameter, private :: ap_val = 0.9_dp ! GeV**-2
  real(dp), parameter, private :: a0_sea = 1.1_dp
  real(dp), parameter, private :: ap_sea = 0.15_dp ! GeV**-2

  ! Power series coefficients; cf. Eq. (16)
  real(dp), parameter, private, dimension(0:4) :: Cu = &
      & [ 1.52_dp, 2.88_dp, -0.095_dp, 0.0_dp, 0.0_dp ]
  real(dp), parameter, private, dimension(0:4) :: Cd = &
      & [ 0.76_dp, 3.11_dp, -3.99_dp, 0.0_dp, 0.0_dp ]
  real(dp), parameter, private, dimension(0:4) :: Cg = &
      & [ 2.23_dp, 5.43_dp, -34.0_dp, 40.6_dp, 0.0_dp ]
  real(dp), parameter, private, dimension(0:4) :: Cs = &
      & [ 0.123_dp, -0.327_dp, 0.692_dp, -0.486_dp, 0.0_dp ]

  ! Misc
  real(dp), parameter, private :: kappa_s = 1.68_dp ! Eq. (20)

  ! ad hoc multiplicative factors to ensure normalizations are exact,
  ! to an absurd number of decimal points
  real(dp), parameter, private :: Nu = 1.0089593_dp
  real(dp), parameter, private :: Nd = 0.9964662_dp
  real(dp), parameter, private :: Ng = 1.0298697_dp

  ! Htilde-type GPDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! from Table 2 of Kroll:2012sm
  ! using Q = Q0, since the code does evolution

  ! Regge
  real(dp), parameter, private :: a0_til = 0.48_dp
  real(dp), parameter, private :: ap_til = 0.45_dp ! GeV**-2

  ! Power series coefficients
  ! Here I assume the form of Eq. (16), but Eq. (28) is actually used.
  ! These are equivalent if odd powers are skipped in Eq. (16).
  real(dp), parameter, private, dimension(0:4) :: Cu_tilde = &
      & [ 0.17_dp, 0.0_dp, 1.34_dp, 0.0_dp, 0.12_dp ]
  real(dp), parameter, private, dimension(0:4) :: Cd_tilde = &
      & [ -0.32_dp, 0.0_dp, -1.427_dp, 0.0_dp, 0.692_dp ]

  real(dp), parameter, private :: eta_u =  0.926_dp ! Eq. (30)
  real(dp), parameter, private :: eta_d = -0.341_dp ! Eq. (30)

  real(dp), parameter, private :: Au =  3.563_dp ! via Eq. (31)
  real(dp), parameter, private :: Ad = -2.528_dp ! via Eq. (31)

  ! Public routines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public :: Hu, Hd, Hs, Hg, Hu_val, Hd_val, Hu_sea, Hd_sea, &
      & Hu_tilde, Hd_tilde

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Direct interfaces: H-type GPDs

    function Hu(x,xi,t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Hu_val(x,xi,t) + Hu_sea(x,xi,t)
    end function Hu

    function Hd(x,xi,t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Hd_val(x,xi,t) + Hd_sea(x,xi,t)
    end function Hd

    function Hs(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = H_n2(x, xi, t, Cs, a0_sea, ap_sea, -1.0_dp)
    end function Hs

    function Hg(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Ng * H_n2(x, xi, t, Cg, a0_sea-1.0_dp, ap_sea, 1.0_dp)
    end function Hg

    function Hu_val(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Nu * H_n1(x, xi, t, Cu, a0_val, ap_val)
    end function Hu_val

    function Hd_val(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Nd * H_n1(x, xi, t, Cd, a0_val, ap_val)
    end function Hd_val

    function Hu_sea(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = kappa_s * H_n2(x, xi, t, Cs, a0_sea, ap_sea, -1.0_dp)
    end function Hu_sea

    function Hd_sea(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = kappa_s * H_n2(x, xi, t, Cs, a0_sea, ap_sea, -1.0_dp)
    end function Hd_sea

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Direct interfaces: Htilde-type GPDs

    function Hu_tilde(x,xi,t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = eta_u*Au*H_n1(x, xi, t, Cu_tilde, a0_til, ap_til)
    end function Hu_tilde

    function Hd_tilde(x,xi,t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = eta_d*Ad*H_n1(x, xi, t, Cd_tilde, a0_til, ap_til)
    end function Hd_tilde

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Switchers

    function H_n1(x, xi, t, C, delta, ap) result(H)
        real(dp), intent(in) :: x, xi, t, delta, ap, C(0:4)
        real(dp) :: H
        !
        if(xi < 0.01*abs(x)) then
          H = H_n1_zero(x, t, C, delta, ap)
        elseif( abs(x - xi) < eps) then
          H = H_n1_analytic(xi+eps, xi, t, C, delta, ap)
        elseif( abs(x + xi) < eps) then
          H = H_n1_analytic(-xi-eps, xi, t, C, delta, ap)
        else
          H = H_n1_analytic(x, xi, t, C, delta, ap)
        endif
    end function H_n1

    function H_n2(x, xi, t, C, delta, ap, sgn) result(H)
        real(dp), intent(in) :: x, xi, t, delta, ap, C(0:4), sgn
        real(dp) :: H
        !
        if(xi < 0.01*abs(x)) then
          H = H_n2_zero(x, t, C, delta, ap, sgn)
        elseif( abs(x - xi) < eps) then
          H = H_n2_analytic(xi+eps, xi, t, C, delta, ap, sgn)
        elseif( abs(x + xi) < eps) then
          H = H_n2_analytic(-xi-eps, xi, t, C, delta, ap, sgn)
        else
          H = H_n2_analytic(x, xi, t, C, delta, ap, sgn)
        endif
    end function H_n2

end module gk

! kernels_nlo.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created December, 2024
!
! All kernels are taken from:
!   Belitsky, Fruend and Mueller
!   Nuclear Physics B 574 (2000) 347-406
!   Belitsky:1999hf
!   arxiv:hep-ph/9912379
!
! About this module:
!
! 1.) TODO
!
! 2.) Stuff independent of the number of active flavors (nfl) and linear in nfl
! have been separated. (TODO)

module kernels_nlo
  use constants,      only: CF, CA, TF, pi, zeta2
  use kernels_common

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  real(dp), parameter, private :: eps = 1e-9_dp

  public :: KV1_Gq_reg, KV1_Gq_reg_nfl
  ! TODO : the rest

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! qq kernel pieces

    ! TODO

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! qG kernel pieces

    ! TODO

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Gq kernel pieces

    function KV1_Gq_reg(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_Gq_half(X1,Y1) - KV1_Gq_half(X2,Y2)
        K = 0.5*K/xi
        K = K * (2.*xi) ! Eq. (15) of BFM
    end function KV1_Gq_reg

    function KV1_Gq_reg_nfl(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_Gq_half_nfl(X1,Y1) - KV1_Gq_half_nfl(X2,Y2)
        K = 0.5*K/xi
        K = K * (2.*xi) ! Eq. (15) of BFM
    end function KV1_Gq_reg_nfl

    function KV1_Gq_half(X, Y) result(K)
        ! Eq. (184)
        real(dp), intent(in) :: X, Y
        real(dp) :: K
        !
        real(dp) :: piece1, piece2, piece3, piece4, piece5
        real(dp) :: GQfV, GQfA, GQfbarV, GQfbarA, GQhV, GQhbarV
        real(dp) :: Xbar, Ybar, beta0
        !real(dp) :: log2X, log2Xbar, logXbarlogX
        K = 0.0_dp
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        ! Some conditions to avoid nans and instability
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        ! The order 1 part of beta0
        beta0 = -11./3.*CA
        ! Terms that contribute
        GQfV    = GQ_f_a(X,Y) + 2.*GQ_f_c(X,Y)
        GQfA    = GQ_f_a(X,Y)
        GQfbarV = GQ_f_a(Xbar,Ybar) + 2.*GQ_f_c(Xbar,Ybar)
        GQfbarA = GQ_f_a(Xbar,Ybar)
        GQhV    = GQ_h_V(X,Y)
        GQhbarV = GQ_hbar_V(Xbar,Y)
        !GQhbarV = GQ_hbar_V(X,Y)
        ! pieces 1-3 have support in both regions
        piece1 = CF**2*( &
            & 3.*GQfA - 4.5*X/Y - (0.5*X*(2.+X)/Ybar - X*Xbar/(Y*Ybar))*abslog(X/Y) &
            & - (1.5*GQfV + 1.5*GQfbarV + X*Xbar/(Y*Ybar))*abslog(1.-X/Y) &
            & - 0.5*GQfA*log2(X/Y) &
            & - 0.5*(GQfV + GQfbarV)*log2(1.-X/Y) &
            & )
        piece2 = -0.5*CF*beta0*( &
            10./3.*GQfV + 2.*X*Xbar/Y + (GQfV+GQfbarV)*abslog(1.-X/Y) &
            & - GQfA*abslog(x) + GQfbarA*abslog(Xbar) &
            & )
        piece3 = CF*CA*( &
            (4./9.-2.*zeta2)*GQfV - 26./9.*X**2/Y + X*(3.+4.*X**2)/Y &
            & - X/9.*(105. - 246.*X + 188.*X**2) &
            & - (X*Xbar/(Y*Ybar) - X*(6.-11.*X+8.*X**2)/Ybar + 4.*X*Xbar*(Xbar-X))*abslog(X/Y) &
            & + X*Xbar/(Y*Ybar)*abslog(1.-X/Y) &
            & + 0.5*(2.*GQfV + 2.*GQfA + GQfbarA - 3.*X**2/Ybar - 2.)*log2(X/Y) &
            & + 0.5*(2. - GQfA - GQfbarA)*log2(Y/X-1.) &
            & - 0.5*GQhV - 0.5*GQhbarV &
            & )
        ! piece 4 and piece 5 have support in ERBL region only
        piece4 = CF**2*( &
            & - 0.5*(X*(2.-5.*X)/Y - 3.*X**2/Ybar)*abslog(X) &
            & + 0.5*GQfA*log2(X) &
            & )
        piece5 = CF*CA*( &
            & (X*Xbar/(Y*Ybar) - X*(6.-11.*X+8.*X**2)/Ybar + X*(18.-63.*X+62.*X**2)/3.)*abslog(X) &
            & - GQfV*logprod(X,Xbar) - X**2*(1.-4.*Y)/(2.*Y*Ybar)*log2(X) &
            & )
        ! Factor in the support regions
        piece1 = piece1*rho_step(X,Y) !* 0.0_dp ! test
        piece2 = piece2*rho_step(X,Y) !* 0.0_dp ! test
        piece3 = piece3*rho_step(X,Y) !* 0.0_dp ! test
        piece4 = piece4*erbl_step(X)  !* 0.0_dp ! test
        piece5 = piece5*erbl_step(X)  !* 0.0_dp ! test
        K = piece1 + piece2 + piece3 + piece4 + piece5
    end function KV1_Gq_half

    function KV1_Gq_half_nfl(X, Y) result(K)
        real(dp), intent(in) :: X, Y
        real(dp) :: K
        !
        real(dp) :: GQfV, GQfA, GQfbarV, GQfbarA
        real(dp) :: Xbar, Ybar, beta0
        K = 0.0_dp
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        ! Some conditions to avoid nans and instability
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        ! The order nfl part of beta0
        beta0 = 4./3.*TF
        ! Terms that contribute
        GQfV    = GQ_f_a(X,Y) + 2.*GQ_f_c(X,Y)
        GQfA    = GQ_f_a(X,Y)
        GQfbarV = GQ_f_a(Xbar,Ybar) + 2.*GQ_f_c(Xbar,Ybar)
        GQfbarA = GQ_f_a(Xbar,Ybar)
        ! Formula
        K = -0.5*CF*beta0*( &
            10./3.*GQfV + 2.*X*Xbar/Y + (GQfV+GQfbarV)*abslog(1.-X/Y) &
            & - GQfA*abslog(x) + GQfbarA*abslog(Xbar) &
            & )
        ! Factor in the support regions
        K = K*rho_step(X,Y)
        !K = 0.0_dp ! test
    end function KV1_Gq_half_nfl

    ! TODO

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GG kernel pieces

    ! TODO

end module kernels_nlo

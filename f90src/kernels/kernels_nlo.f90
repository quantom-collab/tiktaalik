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
  use constants,      only: CF, CA, TF, pi, zeta2, zeta3
  use integration,    only: integrate, integrate2
  use kernels_common

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  real(dp), parameter, private :: eps = 1e-9_dp

  public :: KV1_NSp_pls, KV1_NSp_cst, KV1_NSp_pls_nfl, KV1_NSp_cst_nfl, &
      & KV1_NSm_pls, KV1_NSm_cst, KV1_NSm_pls_nfl, KV1_NSm_cst_nfl, &
      & KV1_qq_reg, &
      & KV1_qG_reg, KV1_Gq_reg, KV1_Gq_reg_nfl, &
      KV1_GG_cst, KV1_GG_pls, KV1_GG_cst_nfl, KV1_GG_pls_nfl
  ! TODO : the rest

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! qq kernel pieeces: nfl-independent terms

    function KV1_NSp_pls(x, y, xi) result(K)
        ! Eq. (177)
        ! Plus-type, plus distribution piece
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_qq_half(X1,Y1,X2,Y2,1.0_dp) + KV1_qq_half(X2,Y2,X1,Y1,1.0_dp)
        K = 0.5*K/xi
    end function KV1_NSp_pls

    function KV1_NSp_cst(x, xi) result(K)
        ! Eq. (177)
        ! Plus-type, constant term
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        ! Explicit term in Eq. (177)
        K = CF*(CF - 0.5*CA)*(6.5-6.*zeta2+4.*zeta3)
        ! From plus prescription
        K = K + integrate2(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = KV1_NSp_pls(x,y,xi) - KV1_NSp_pls(y,x,xi)
          end function integrand
    end function KV1_NSp_cst

    function KV1_NSm_pls(x, y, xi) result(K)
        ! Eq. (177)
        ! Minus-type, plus distribution piece
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_qq_half(X1,Y1,X2,Y2,-1.0_dp) + KV1_qq_half(X2,Y2,X1,Y1,-1.0_dp)
        K = 0.5*K/xi
    end function KV1_NSm_pls

    function KV1_NSm_cst(x, xi) result(K)
        ! Eq. (177)
        ! Minus-type, constant term
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        ! Constant only from plus prescription
        K = integrate2(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = KV1_NSm_pls(x,y,xi) - KV1_NSm_pls(y,x,xi)
          end function integrand
    end function KV1_NSm_cst

    function KV1_qq_half(X, Y, Xbar, Ybar, zsign) result(K)
        ! Eq. (177)
        real(dp), intent(in) :: X, Y, Xbar, Ybar, zsign
        real(dp) :: K
        !
        real(dp) :: QQf, QQfbar, QQh, QQhbar, beta0
        real(dp) :: piece1, piece2, piece3, piece4
        K = 0.0_dp
        ! Some conditions to avoid nans
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        ! Auxiliary quantities
        QQf    = QQ_f_a(X   ,Y)    + QQ_f_b(X   ,Y   )
        QQfbar = QQ_f_a(Xbar,Ybar) + QQ_f_b(Xbar,Ybar)
        QQh = QQ_h(X,Y)
        QQhbar = QQ_hbar(Xbar,Y)
        beta0 = -11./3.*CA
        ! pieces 1-3 have support in both the ERBL and DLGAP regions
        piece1 = CF**2*( (4./3.-2.*zeta2)*QQf + 3.*X/Y &
            & - (1.5*QQf-0.5*X/Ybar)*abslog(X/Y) &
            & - (QQf-QQfbar)*logprod(X/Y,1.-X/Y) &
            & + (QQf + 0.5*X/Ybar)*log2(X/Y) )
        piece2 = -0.5*CF*beta0*( 5./3.*QQf + X/Y + QQf*abslog(X/Y) )
        piece3 = - CF*(CF - 0.5*CA)*( 4./3.*QQf + 2.*X/Y + QQh - zsign*QQhbar )
        ! piece 4 only has support in the ERBL region
        piece4 = -0.5*X/Ybar*(abslog(X) + log2(X) - 2.*logprod(X,Xbar))*CF**2
        ! Factor in the support regions
        piece1 = piece1*rho_step(X,Y)
        piece2 = piece2*rho_step(X,Y)
        piece3 = piece3*rho_step(X,Y)
        piece4 = piece4*erbl_step(X)
        K = piece1 + piece2 + piece3 + piece4
    end function KV1_qq_half

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! qq kernel pieces: terms linear in nfl

    function KV1_NSp_pls_nfl(x, y, xi) result(K)
        ! Eq. (177)
        ! Plus-type, plus distribution piece
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_qq_half_nfl(X1,Y1,X2,Y2,1.0_dp) + KV1_qq_half_nfl(X2,Y2,X1,Y1,1.0_dp)
        K = 0.5*K/xi
    end function KV1_NSp_pls_nfl

    function KV1_NSp_cst_nfl(x, xi) result(K)
        ! Eq. (177)
        ! Plus-type, constant term
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        ! From plus prescription
        K = integrate2(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = KV1_NSp_pls_nfl(x,y,xi) - KV1_NSp_pls_nfl(y,x,xi)
          end function integrand
    end function KV1_NSp_cst_nfl

    function KV1_NSm_pls_nfl(x, y, xi) result(K)
        ! Eq. (177)
        ! Minus-type, plus distribution piece
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_qq_half_nfl(X1,Y1,X2,Y2,-1.0_dp) + KV1_qq_half_nfl(X2,Y2,X1,Y1,-1.0_dp)
        K = 0.5*K/xi
    end function KV1_NSm_pls_nfl

    function KV1_NSm_cst_nfl(x, xi) result(K)
        ! Eq. (177)
        ! Minus-type, constant term
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        ! Constant only from plus prescription
        K = integrate2(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = KV1_NSm_pls_nfl(x,y,xi) - KV1_NSm_pls_nfl(y,x,xi)
          end function integrand
    end function KV1_NSm_cst_nfl

    function KV1_qq_half_nfl(X, Y, Xbar, Ybar, zsign) result(K)
        ! Eq. (177)
        real(dp), intent(in) :: X, Y, Xbar, Ybar, zsign
        real(dp) :: K
        !
        real(dp) :: QQf, QQfbar, QQh, QQhbar, beta0
        !!!real(dp) :: piece1, piece2, piece3, piece4
        K = 0.0_dp
        ! Some conditions to avoid nans
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        ! Auxiliary quantities
        QQf    = QQ_f_a(X   ,Y)    + QQ_f_b(X   ,Y   )
        QQfbar = QQ_f_a(Xbar,Ybar) + QQ_f_b(Xbar,Ybar)
        QQh = QQ_h(X,Y)
        QQhbar = QQ_hbar(Xbar,Y)
        beta0 = 4./3.*TF
        K = -0.5*CF*beta0*( 5./3.*QQf + X/Y + QQf*abslog(X/Y) )
        K = K*rho_step(X,Y)
    end function KV1_qq_half_nfl

    function KV1_qq_reg(x, y, xi) result(K)
        ! Note: the singelt QQ kernel is the plus-type NS kernel plus this!
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_qq_half_sin(X1,Y1) + KV1_qq_half_sin(X2,Y2)
        K = 0.5*K/xi
    end function KV1_qq_reg

    function KV1_qq_half_sin(X, Y) result(K)
        ! Entire thing needs to be multiplied by nfl
        real(dp), intent(in) :: X, Y
        real(dp) :: K
        !
        real(dp) :: piece1, piece2
        real(dp) :: Xbar, Ybar
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        ! piece 1 has support in both ERBL and DGLAP regions
        piece1 = X*(3.-8.*X*Ybar)/Y + X*(5.-8.*X)/Ybar*abslog(X/Y) &
            & + (X/Ybar - 4.*X*Xbar)*log2(X/Y) &
            & + 8.*X*Xbar*(li2(Xbar) - li2(Ybar) + logprod(Xbar,Y))
        ! piece 2 only has support in the ERBL region
        piece2 = -290./9.*X*Xbar &
            & - (X*(5.-8.*X)/Ybar + 2.*X*(9.-19.*Xbar)/3.)*abslog(X) &
            & - (X/Ybar - 4.*X*Xbar)*log2(X) &
            & + 4.*X*Xbar*logprod(X,Xbar)
        ! Factor in the support regions
        piece1 = piece1*rho_step(X,Y)
        piece2 = piece2*erbl_step(X)
        K = 2.*CF*TF*(piece1 + piece2)
    end function KV1_qq_half_sin

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! qG kernel pieces

    function KV1_qG_reg(x, y, xi) result(K)
        ! NOTE: the entire thing should be multiplied by nfl
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_qG_half(X1,Y1) - KV1_qG_half(X2,Y2)
        K = 0.5*K/xi
    end function KV1_qG_reg

    function KV1_qG_half(X, Y) result(K)
        ! Eq. (183)
        real(dp), intent(in) :: X, Y
        real(dp) :: K
        !
        real(dp) :: piece1, piece2, piece3, piece4
        real(dp) :: QGfV, QGfA, QGfbarV, QGfbarA, QGhV, QGhbarV
        real(dp) :: Xbar, Ybar
        K = 0.0_dp
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        ! Some conditions to avoid nans and instability
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        ! Terms that contribute
        QGfV    = -QG_f_a(X,Y) - 2.*QG_f_c(X,Y)
        QGfA    = -QG_f_a(X,Y)
        QGfbarV = -QG_f_a(Xbar,Ybar) - 2.*QG_f_c(Xbar,Ybar)
        QGfbarA = -QG_f_a(Xbar,Ybar)
        QGhV    = QG_h_V(X,Y)
        QGhbarV = QG_hbar_V(Xbar,Y)
        ! piece 1 and piece 2 have support in both regions
        piece1 = TF*CF*( &
            & 2.*(5.-2.*zeta2)*QGfV - 4.*QGfA + X/(Y*Ybar) &
            & + 2.*(QGfV + QGfbarV + (2.*X-2.*Y+X*Y)/(2.*Y*Ybar**2))*abslog(X/Y) &
            & - 2.*(QGfV + QGfbarV + 1./(Y*Ybar))*abslog(1.-X/Y) &
            & + QGfbarA*log2(X/Y) + (QGfV+QGfbarV)*log2(Y/X-1.) &
            & )
        piece2 = TF*CA*( &
            & 6.*X*(3.-4.*X)/Ybar - 2.*X*(11.-16.*X)/Y + 4.*X*(1.-3.*X)/Y**2 &
            & - 2.*(QGfV - 5.*QGfbarV + (5.-7.*X)/Ybar**2 + 2.*X*(3.-2.*X-12.*Xbar*Y)/(Y*Ybar))*abslog(X/Y) &
            & + 2.*(QGfV + QGfbarV + 1./(Y*Ybar))*abslog(1.-X/Y) &
            & + (QGfV - QGfbarV + (1.-4.*X)/Ybar**2)*log2(X/Y) &
            & - (QGfV + QGfbarV)*abslog(1.-X/Y)**2 - QGhV + QGhbarV &
            & )
        ! piece 3 and piece 4 have support in ERBL region only
        piece3 = TF*CF*( &
            & - 2.*(QGfV + QGfbarV + QGfbarA + X*(Y+4.*Ybar**2)/(2.*Y*Ybar**2))*abslog(X) &
            & - 2.*qGfV*logprod(Xbar,X) &
            & + (QGfV - QGfbarV - QGfbarA)*log2(X) &
            & )
        piece4 = TF*CA*( &
            & 2.*(QGfV - 5.*QGfbarV + (5.-7.*X)/Ybar**2 + X*(5.-6.*X-22.*Y+28.*X*Y)/(Y*Ybar))*abslog(X) &
            & - (QGfV - QGfbarV + (1.-4.*X)/Ybar**2)*log2(X) &
            & )
        ! Factor in the support regions
        piece1 = piece1*rho_step(X,Y)
        piece2 = piece2*rho_step(X,Y)
        piece3 = piece3*erbl_step(X)
        piece4 = piece4*erbl_step(X)
        K = piece1 + piece2 + piece3 + piece4
    end function KV1_qG_half

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Gq kernel pieces: nfl-independent

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
        piece1 = piece1*rho_step(X,Y)
        piece2 = piece2*rho_step(X,Y)
        piece3 = piece3*rho_step(X,Y)
        piece4 = piece4*erbl_step(X)
        piece5 = piece5*erbl_step(X)
        K = piece1 + piece2 + piece3 + piece4 + piece5
    end function KV1_Gq_half

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Gq kernel pieces: linear in nfl

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

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GG kernel pieces: nfl-independent

    function KV1_GG_pls(x, y, xi) result(K)
        ! Eq. (177)
        ! Plus-type, plus distribution piece
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_GG_half(X1,Y1) + KV1_GG_half(X2,Y2)
        K = 0.5*K/xi
    end function KV1_GG_pls

    function KV1_GG_cst(x, xi) result(K)
        ! Eq. (177)
        ! Plus-type, constant term
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        ! From plus prescription
        K = 0.0_dp
        K = integrate(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = KV1_GG_pls(x,y,xi) - KV1_GG_pls(y,x,xi)
          end function integrand
    end function KV1_GG_cst

    function KV1_GG_half(X, Y) result(K)
        ! Eq. (185)
        real(dp), intent(in) :: X, Y
        real(dp) :: K
        !
        real(dp) :: piece1, piece2, piece3, piece4, piece5
        real(dp) :: Xbar, Ybar, beta0
        real(dp) :: GGfV, GGfbarV, GGfc, GGhV, GGhbarV
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        K = 0.0_dp
        ! Some conditions to avoid nans and instability
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        beta0 = -11./3.*CA
        ! Terms that contribute
        GGfV    = 2.*GG_f_a(X   ,Y   ) + GG_f_b(X   ,Y   ) + 2.*GG_f_C(X   ,Y   )
        GGfbarV = 2.*GG_f_a(Xbar,Ybar) + GG_f_b(Xbar,Ybar) + 2.*GG_f_C(Xbar,Ybar)
        GGfc    = GG_f_c(X,Y)
        GGhV    = GG_h_V(X,Y)
        GGhbarV = GG_hbar_V(Xbar,Y)
        ! pieces 1, 2, 3 have support in both regions
        piece1 = CA**2*( &
            & (2./3.-2.*zeta2)*GGfV - 0.5*GGfc + X*(4.+5.*X)/(4.*Y**2) &
            & - X*(10.-19.*X+16.*X**2)/Y + X*(6.-13.*X+8.*X**2)/Ybar &
            & + (2.*X*(1.-2.*Xbar*(X*Y+Xbar*Ybar))/(Y*Ybar) + X*(2.-3.*X-8.*Xbar**2)/Ybar**2)*abslog(X/Y) &
            & + (GGfV + 2.*X**2/Ybar**2)*log2(X/Y) &
            & - (GGfV-GGfbarV)*logprod(X/Y,1.-X/Y) &
            & -0.5*GGhV - 0.5*GGhbarV &
            & )
        piece2 = 0.5*CA*beta0*( &
            & - 5./3.*GGfV - 13./3.*GGfc - 11./2.*X**2/Y**2 + X*(Xbar+Ybar)/(Y**2*Ybar) &
            & + X**2/Ybar**2*abslog(Y) + Xbar**2/Y**2*abslog(Xbar) &
            & )
        ! pieces 4, 5 have support only in the ERBL region
        piece4 = CA**2*( &
            & - 2.*X**2/Ybar**2*log2(X) &
            & + (Xbar*(X**2+Xbar**2-2.*Xbar*(1.+2.*X)*Ybar) + Ybar)/Ybar**2*logprod(X,Xbar) &
            & + ((3.*X*(X**2+Xbar**2)*(Ybar-Y)-4.*X**2*Y)/(Y*Ybar) - X*(2.-3.*X-8.*Xbar**2)/Ybar**2)*abslog(X) &
            & )
        ! Factor in the support regions
        piece1 = piece1*rho_step(X,Y)
        piece2 = piece2*rho_step(X,Y)
        piece3 = 0.0_dp ! linear in nfl
        piece4 = piece4*erbl_step(X)
        piece5 = 0.0_dp ! linear in nfl
        K = piece1 + piece2 + piece3 + piece4 + piece5
    end function KV1_GG_half

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GG kernel pieces: linear in nfl

    function KV1_GG_pls_nfl(x, y, xi) result(K)
        ! Eq. (177)
        ! Plus-type, plus distribution piece
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        real(dp) :: X1, X2, Y1, Y2
        X1 = 0.5*(1.+x/xi)
        X2 = 0.5*(1.-x/xi)
        Y1 = 0.5*(1.+y/xi)
        Y2 = 0.5*(1.-y/xi)
        K = KV1_GG_half_nfl(X1,Y1) + KV1_GG_half_nfl(X2,Y2)
        K = 0.5*K/xi
    end function KV1_GG_pls_nfl

    function KV1_GG_cst_nfl(x, xi) result(K)
        ! Eq. (177)
        ! Plus-type, constant term
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        ! Explicit term in Eqs. (175) and (176)
        K = -1./108.*(35.*CA + 74.*CF)
        ! From plus prescription
        ! test...
        K = K + integrate(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = KV1_GG_pls_nfl(x,y,xi) - KV1_GG_pls_nfl(y,x,xi)
          end function integrand
    end function KV1_GG_cst_nfl

    function KV1_GG_half_nfl(X, Y) result(K)
        real(dp), intent(in) :: X, Y
        real(dp) :: K
        !
        real(dp) :: piece1, piece2, piece3, piece4, piece5
        real(dp) :: Xbar, Ybar, beta0
        real(dp) :: GGfV, GGfbarV, GGfc!, GGhV, GGhbarV
        Xbar = 1.0_dp - X
        Ybar = 1.0_dp - Y
        K = 0.0_dp
        ! Some conditions to avoid nans and instability
        if(abs(X-Y)  < eps) return
        if(abs(Y)    < eps) return
        if(abs(Ybar) < eps) return
        if(abs(X)    < eps) return
        if(abs(Xbar) < eps) return
        beta0 = 4./3.*TF
        ! Terms that contribute
        GGfV    = 2.*GG_f_a(X   ,Y   ) + GG_f_b(X   ,Y   ) + 2.*GG_f_C(X   ,Y   )
        GGfbarV = 2.*GG_f_a(Xbar,Ybar) + GG_f_b(Xbar,Ybar) + 2.*GG_f_C(Xbar,Ybar)
        GGfc    = GG_f_c(X,Y)
        !GGhV    = GG_h_V(X,Y)
        !GGhbarV = GG_hbar_V(Xbar,Y)
        ! pieces 1, 2, 3 have support in both regions
        piece2 = 0.5*CA*beta0*( &
            & - 5./3.*GGfV - 13./3.*GGfc - 11./2.*X**2/Y**2 + X*(Xbar+Ybar)/(Y**2*Ybar) &
            & + X**2/Ybar**2*abslog(Y) + Xbar**2/Y**2*abslog(Xbar) &
            & )
        piece3 = CF*TF*( &
            & -20./3.*GGfc - 12.*X**2/Y**2 + 2.*X*(2.+Y+3.*X*Y-4.*Y**2)/(Y**2*Ybar) &
            & - 2.*X*(X-Y+2.*X*Y)/(Y*Ybar**2)*abslog(X/Y) - X**2/Ybar**2*log2(X/Y) &
            & )
        ! pieces 4, 5 have support only in the ERBL region
        piece5 = CF*TF*( &
            & -2.*X*(1.-3.*X-2.*X*Ybar)/Ybar**2*abslog(X) + X**2/Ybar**2*log2(X) &
            & )
        ! Factor in the support regions
        piece1 = 0.0_dp ! nfl-independent
        piece2 = piece2*rho_step(X,Y)
        piece3 = piece3*rho_step(X,Y)
        piece4 = 0.0_dp ! nfl-independent
        piece5 = piece5*erbl_step(X)
        K = piece1 + piece2 + piece3 + piece4 + piece5
    end function KV1_GG_half_nfl

end module kernels_nlo

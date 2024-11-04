module integration

  implicit none
  private

  integer,  parameter :: dp = kind(1d0)

  public :: integrate

  contains

    ! ==========================================================================
    ! Integration method used by tiktaalik
    ! Uses 15-point Gauss-Kronrod rule in six different regions of y domain.
    ! Calls qk15 from quadpack.

    function integrate(func, x, xi) result(integral)
        real(dp), external :: func
        real(dp), intent(in) :: x, xi
        real(dp) :: integral
        !
        real(dp), parameter :: eps = 1e-12_dp ! to avoid undefined behavior
        real(dp) :: aa(7), ii(6), abserr(6), resabs(6), resasc(6)
        integer :: i
        aa(1) = -1.0_dp
        aa(2) = -abs(x)
        aa(3) = -abs(xi)
        aa(4) = 0.0_dp
        aa(5) = abs(xi)
        aa(6) = abs(x)
        aa(7) = 1.0_dp
        if(aa(3) < aa(2)) call swap(aa(2),aa(3))
        if(aa(6) < aa(5)) call swap(aa(5),aa(6))
        do i=1, 6, 1
          call qk15(func, aa(i)+eps, aa(i+1)-eps, ii(i), abserr(i), resabs(i), resasc(i))
        end do
        integral = sum(ii)
    end function integrate

    subroutine swap(x,y)
        real(dp), intent(inout) :: x, y
        real(dp) :: z
        z = x
        x = y
        y = z
    end subroutine swap

    ! ==========================================================================
    ! qk15 routine from quadpack

    subroutine qk15(f,a,b,result,abserr,resabs,resasc)
        !! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
        !
        ! This was originally from the quadpack package.
        ! Only this subrotuine is used by tiktaalik, so this subroutine has been isolated.
        ! Change by Ian Cloet: dp parameter used to control float precision.
        !
        ! Discussion:
        !
        !   This routine approximates
        !     I = integral ( A <= X <= B ) F(X) dx
        !   with an error estimate, and
        !     J = integral ( A <= X <= B ) | F(X) | dx
        !
        ! Reference:
        !
        !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
        !   QUADPACK, a Subroutine Package for Automatic Integration,
        !   Springer Verlag, 1983
        !
        ! Parameters:
        !
        !   Input, external real F, the name of the function routine, of the form
        !     function f ( x )
        !     real f
        !     real x
        !   which evaluates the integrand function.
        !
        !   Input, real A, B, the limits of integration.
        !
        !   Output, real RESULT, the estimated value of the integral.
        !   RESULT is computed by applying the 15-point Kronrod rule (RESK)
        !   obtained by optimal addition of abscissae to the 7-point Gauss rule
        !   (RESG).
        !
        !   Output, real ABSERR, an estimate of | I - RESULT |.
        !
        !   Output, real RESABS, approximation to the integral of the absolute
        !   value of F.
        !
        !   Output, real RESASC, approximation to the integral | F-I/(B-A) |
        !   over [A,B].
        !
        ! Local Parameters:
        !
        !          the abscissae and weights are given for the interval (-1,1).
        !          because of symmetry only the positive abscissae and their
        !          corresponding weights are given.
        !
        !          xgk    - abscissae of the 15-point Kronrod rule
        !                   xgk(2), xgk(4), ...  abscissae of the 7-point
        !                   Gauss rule
        !                   xgk(1), xgk(3), ...  abscissae which are optimally
        !                   added to the 7-point Gauss rule
        !
        !          wgk    - weights of the 15-point Kronrod rule
        !
        !          wg     - weights of the 7-point Gauss rule
        !
        !          centr  - mid point of the interval
        !          hlgth  - half-length of the interval
        !          absc   - abscissa
        !          fval*  - function value
        !          resg   - result of the 7-point Gauss formula
        !          resk   - result of the 15-point Kronrod formula
        !          reskh  - approximation to the mean value of f over (a,b),
        !                   i.e. to i/(b-a)

        real(dp), external    :: f
        real(dp), intent(in)  :: a,b
        real(dp), intent(out) :: result,abserr,resabs,resasc

        ! Declare local variables
        integer  :: j,jtw,jtwm1
        real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(7),fv2(7),hlgth, &
                    resg,resk,reskh,wg(4),wgk(8),xgk(8)

        ! TODO: these can probably be made constants at compile time
        xgk = (/ 9.914553711208126e-01_dp,   9.491079123427585e-01_dp, &
                 8.648644233597691e-01_dp,   7.415311855993944e-01_dp, &
                 5.860872354676911e-01_dp,   4.058451513773972e-01_dp, &
                 2.077849550078985e-01_dp,   0.0_dp                  /)
        wgk = (/ 2.293532201052922e-02_dp,   6.309209262997855e-02_dp, &
                 1.047900103222502e-01_dp,   1.406532597155259e-01_dp, &
                 1.690047266392679e-01_dp,   1.903505780647854e-01_dp, &
                 2.044329400752989e-01_dp,   2.094821410847278e-01_dp    /)
        wg  = (/ 1.294849661688697e-01_dp,   2.797053914892767e-01_dp, &
                 3.818300505051189e-01_dp,   4.179591836734694e-01_dp    /)

        centr  = 5.0e-01_dp*(a+b)
        hlgth  = 5.0e-01_dp*(b-a)
        dhlgth = abs(hlgth)

        ! Compute the 15-point Kronrod approximation to the integral,
        ! and estimate the absolute error.
        fc = f(centr)
        resg = fc*wg(4)
        resk = fc*wgk(8)
        resabs = abs(resk)

        do j = 1, 3
           jtw = j*2
           absc = hlgth*xgk(jtw)
           fval1 = f(centr-absc)
           fval2 = f(centr+absc)
           fv1(jtw) = fval1
           fv2(jtw) = fval2
           fsum = fval1+fval2
           resg = resg+wg(j)*fsum
           resk = resk+wgk(jtw)*fsum
           resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
        enddo

        do j = 1, 4
           jtwm1 = j*2-1
           absc = hlgth*xgk(jtwm1)
           fval1 = f(centr-absc)
           fval2 = f(centr+absc)
           fv1(jtwm1) = fval1
           fv2(jtwm1) = fval2
           fsum = fval1+fval2
           resk = resk+wgk(jtwm1)*fsum
           resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
        enddo

        reskh = resk * 5.0e-01_dp
        resasc = wgk(8)*abs(fc-reskh)

        do j = 1, 7
           resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
        enddo

        result = resk*hlgth
        resabs = resabs*dhlgth
        resasc = resasc*dhlgth
        abserr = abs((resk-resg)*hlgth)

        if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp ) &
           abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

        if ( resabs > tiny ( resabs ) / (5.0e+01_dp* epsilon ( resabs ) ) ) &
           abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

   end subroutine qk15

end module integration

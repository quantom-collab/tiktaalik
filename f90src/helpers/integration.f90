module integration
  use gridspace, only: push_forward, pull_back, push_jacob

  implicit none
  private

  integer,  parameter :: dp = kind(1d0)

  public :: integrate, adaptive_integrate

  contains

    ! ==========================================================================
    ! Integration method used by tiktaalik
    ! Uses 21-point Gauss-Kronrod rule in six different regions of y domain.
    ! Calls qk21 from quadpack.

    function integrate(func, x, xi, n_pixels, grid_type) result(integral)
        ! Do the spacing in eta instead, but break into regions the same way
        real(dp), external   :: func
        real(dp), intent(in) :: x, xi
        integer,  intent(in) :: n_pixels, grid_type
        real(dp) :: integral
        !
        real(dp) :: aa(7), ii(6), abserr(6), resabs(6), resasc(6)
        integer :: i
        aa(1) = pull_back( -1.0_dp, xi, n_pixels, grid_type)
        aa(2) = pull_back( -abs(x), xi, n_pixels, grid_type)
        aa(3) = pull_back(-abs(xi), xi, n_pixels, grid_type)
        aa(4) = pull_back  (0.0_dp, xi, n_pixels, grid_type)
        aa(5) = pull_back( abs(xi), xi, n_pixels, grid_type)
        aa(6) = pull_back(  abs(x), xi, n_pixels, grid_type)
        aa(7) = pull_back(  1.0_dp, xi, n_pixels, grid_type)
        if(aa(3) < aa(2)) call swap(aa(2),aa(3))
        if(aa(6) < aa(5)) call swap(aa(5),aa(6))
        do i=1, 6, 1
          call qk21(pulled_back_func, aa(i), aa(i+1), ii(i), abserr(i), resabs(i), resasc(i))
        end do
        integral = sum(ii)
        return
        contains
          function pulled_back_func(eta) result(f)
              real(dp), intent(in) :: eta
              real(dp) :: f, J, x_
              x_ = push_forward(eta, xi, n_pixels, grid_type)
              J  = push_jacob(  eta, xi, n_pixels, grid_type)
              f  = func(x_) * J
          end function pulled_back_func
    end function integrate

    !function integrate_depr(func, x, xi) result(integral)
    !    real(dp), external :: func
    !    real(dp), intent(in) :: x, xi
    !    real(dp) :: integral
    !    !
    !    real(dp) :: aa(7), ii(6), abserr(6), resabs(6), resasc(6)
    !    integer :: i
    !    aa(1) = -1.0_dp
    !    aa(2) = -abs(x)
    !    aa(3) = -abs(xi)
    !    aa(4) = 0.0_dp
    !    aa(5) = abs(xi)
    !    aa(6) = abs(x)
    !    aa(7) = 1.0_dp
    !    if(aa(3) < aa(2)) call swap(aa(2),aa(3))
    !    if(aa(6) < aa(5)) call swap(aa(5),aa(6))
    !    do i=1, 6, 1
    !      call qk21(func, aa(i), aa(i+1), ii(i), abserr(i), resabs(i), resasc(i))
    !    end do
    !    integral = sum(ii)
    !end function integrate_depr

    function adaptive_integrate(func, x, xi) result(integral)
        real(dp), external :: func
        real(dp), intent(in) :: x, xi
        real(dp) :: integral
        !
        real(dp), parameter :: eps  =  1e-12_dp ! to avoid undefined behavior
        real(dp), parameter :: ymin = -1.0_dp + eps
        real(dp), parameter :: ymax =  1.0_dp - eps
        integral = 0.0_dp
        if(abs(x-xi) < eps) then
          integral = 0.0_dp
        elseif(abs(x+xi) < eps) then
          integral = 0.0_dp
        elseif(x > xi) then
          integral = integral + iqags(func,  ymin,        x-eps)
          integral = integral + iqags(func,  x+eps,       ymax)
        elseif(x < -xi) then
          integral = integral + iqags(func,  ymin,        x-eps)
          integral = integral + iqags(func,  x+eps,       ymax)
        elseif(x > -xi .and. x < xi) then
          integral = integral + iqags(func,  ymin,       -abs(x)-eps)
          integral = integral + iqags(func, -abs(x)+eps,  abs(x)-eps)
          integral = integral + iqags(func,  abs(x)+eps,  ymax)
        endif
    end function adaptive_integrate

    subroutine swap(x,y)
        real(dp), intent(inout) :: x, y
        real(dp) :: z
        z = x
        x = y
        y = z
    end subroutine swap

    ! ==========================================================================
    ! qk21 routine from quadpack

    subroutine qk21(f,a,b,result,abserr,resabs,resasc)
        ! QK21 carries out a 21-point Gauss-Kronrod quadrature rule.
        !
        ! This was originally from the quadpack package.
        ! Change by Ian Cloet: dp parameter used to control float precision.
        !
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
        !                      result is computed by applying the 21-point
        !                      Kronrod rule (resk) obtained by optimal addition
        !                      of abscissae to the 10-point Gauss rule (resg).
        !
        !   Output, real ABSERR, an estimate of | I - RESULT |.
        !
        !   Output, real RESABS, approximation to the integral of the absolute
        !   value of F.
        !
        !   Output, real RESASC, approximation to the integral | F-I/(B-A) |
        !   over [A,B].

        real(dp), external    :: f
        real(dp), intent(in)  :: a,b
        real(dp), intent(out) :: result,abserr,resabs,resasc

        ! Declare local variables
        integer  :: j,jtw,jtwm1
        real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(10),fv2(10),hlgth, &
            & resg,resk,reskh,wg(5),wgk(11),xgk(11)

        !   the abscissae and weights are given for the interval (-1,1).
        !   because of symmetry only the positive abscissae and their
        !   corresponding weights are given.
        !   xgk    - abscissae of the 21-point Kronrod rule
        !            xgk(2), xgk(4), ...  abscissae of the 10-point
        !            Gauss rule
        !            xgk(1), xgk(3), ...  abscissae which are optimally
        !            added to the 10-point Gauss rule
        !   wgk    - weights of the 21-point Kronrod rule
        !   wg     - weights of the 10-point Gauss rule

        ! TODO: these should be made constant at compile time
        xgk = (/ 9.956571630258081e-01_dp,     9.739065285171717e-01_dp, &
               & 9.301574913557082e-01_dp,     8.650633666889845e-01_dp, &
               & 7.808177265864169e-01_dp,     6.794095682990244e-01_dp, &
               & 5.627571346686047e-01_dp,     4.333953941292472e-01_dp, &
               & 2.943928627014602e-01_dp,     1.488743389816312e-01_dp, &
               & 0.000000000000000_dp                               /)
        wgk = (/ 1.169463886737187e-02_dp,     3.255816230796473e-02_dp, &
               & 5.475589657435200e-02_dp,     7.503967481091995e-02_dp, &
               & 9.312545458369761e-02_dp,     1.093871588022976e-01_dp, &
               & 1.234919762620659e-01_dp,     1.347092173114733e-01_dp, &
               & 1.427759385770601e-01_dp,     1.477391049013385e-01_dp, &
               & 1.494455540029169e-01_dp                               /)
        wg  = (/ 6.667134430868814e-02_dp,     1.494513491505806e-01_dp, &
               & 2.190863625159820e-01_dp,     2.692667193099964e-01_dp, &
               & 2.955242247147529e-01_dp                               /)

        !   list of major variables
        !
        !   centr  - mid point of the interval
        !   hlgth  - half-length of the interval
        !   absc   - abscissa
        !   fval*  - function value
        !   resg   - result of the 10-point Gauss formula
        !   resk   - result of the 21-point Kronrod formula
        !   reskh  - approximation to the mean value of f over (a,b),
        !            i.e. to i/(b-a)
        centr  = 5.0e-01_dp*(a+b)
        hlgth  = 5.0e-01_dp*(b-a)
        dhlgth = abs(hlgth)

        ! Compute the 21-point Kronrod approximation to the
        ! integral, and estimate the absolute error.
        resg = 0.0_dp
        fc = f(centr)
        resk = wgk(11)*fc
        resabs = abs(resk)

        do j = 1, 5
          jtw = 2*j
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

        do j = 1, 5
          jtwm1 = 2*j-1
          absc = hlgth*xgk(jtwm1)
          fval1 = f(centr-absc)
          fval2 = f(centr+absc)
          fv1(jtwm1) = fval1
          fv2(jtwm1) = fval2
          fsum = fval1+fval2
          resk = resk+wgk(jtwm1)*fsum
          resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
        enddo

        reskh = resk*5.0e-01_dp
        resasc = wgk(11)*abs(fc-reskh)

        do j = 1, 10
          resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
        enddo

        result = resk*hlgth
        resabs = resabs*dhlgth
        resasc = resasc*dhlgth
        abserr = abs((resk-resg)*hlgth)

        if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp) then
          abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)
        endif

        if ( resabs > tiny ( resabs ) /(5.0e+01_dp* epsilon ( resabs ) )) then
          abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)
        endif

    end subroutine qk21

    ! ==========================================================================
    ! qk15 routine from quadpack

    !subroutine qk15(f,a,b,result,abserr,resabs,resasc)
    !    !! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
    !    !
    !    ! This was originally from the quadpack package.
    !    ! Only this subrotuine is used by tiktaalik, so this subroutine has been isolated.
    !    ! Change by Ian Cloet: dp parameter used to control float precision.
    !    !
    !    ! Discussion:
    !    !
    !    !   This routine approximates
    !    !     I = integral ( A <= X <= B ) F(X) dx
    !    !   with an error estimate, and
    !    !     J = integral ( A <= X <= B ) | F(X) | dx
    !    !
    !    ! Reference:
    !    !
    !    !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
    !    !   QUADPACK, a Subroutine Package for Automatic Integration,
    !    !   Springer Verlag, 1983
    !    !
    !    ! Parameters:
    !    !
    !    !   Input, external real F, the name of the function routine, of the form
    !    !     function f ( x )
    !    !     real f
    !    !     real x
    !    !   which evaluates the integrand function.
    !    !
    !    !   Input, real A, B, the limits of integration.
    !    !
    !    !   Output, real RESULT, the estimated value of the integral.
    !    !   RESULT is computed by applying the 15-point Kronrod rule (RESK)
    !    !   obtained by optimal addition of abscissae to the 7-point Gauss rule
    !    !   (RESG).
    !    !
    !    !   Output, real ABSERR, an estimate of | I - RESULT |.
    !    !
    !    !   Output, real RESABS, approximation to the integral of the absolute
    !    !   value of F.
    !    !
    !    !   Output, real RESASC, approximation to the integral | F-I/(B-A) |
    !    !   over [A,B].
    !    !
    !    ! Local Parameters:
    !    !
    !    !          the abscissae and weights are given for the interval (-1,1).
    !    !          because of symmetry only the positive abscissae and their
    !    !          corresponding weights are given.
    !    !
    !    !          xgk    - abscissae of the 15-point Kronrod rule
    !    !                   xgk(2), xgk(4), ...  abscissae of the 7-point
    !    !                   Gauss rule
    !    !                   xgk(1), xgk(3), ...  abscissae which are optimally
    !    !                   added to the 7-point Gauss rule
    !    !
    !    !          wgk    - weights of the 15-point Kronrod rule
    !    !
    !    !          wg     - weights of the 7-point Gauss rule
    !    !
    !    !          centr  - mid point of the interval
    !    !          hlgth  - half-length of the interval
    !    !          absc   - abscissa
    !    !          fval*  - function value
    !    !          resg   - result of the 7-point Gauss formula
    !    !          resk   - result of the 15-point Kronrod formula
    !    !          reskh  - approximation to the mean value of f over (a,b),
    !    !                   i.e. to i/(b-a)

    !    real(dp), external    :: f
    !    real(dp), intent(in)  :: a,b
    !    real(dp), intent(out) :: result,abserr,resabs,resasc

    !    ! Declare local variables
    !    integer  :: j,jtw,jtwm1
    !    real(dp) :: absc,centr,dhlgth,fc,fsum,fval1,fval2,fv1(7),fv2(7),hlgth, &
    !                resg,resk,reskh,wg(4),wgk(8),xgk(8)

    !    ! TODO: these can probably be made constants at compile time
    !    xgk = (/ 9.914553711208126e-01_dp,   9.491079123427585e-01_dp, &
    !             8.648644233597691e-01_dp,   7.415311855993944e-01_dp, &
    !             5.860872354676911e-01_dp,   4.058451513773972e-01_dp, &
    !             2.077849550078985e-01_dp,   0.0_dp                  /)
    !    wgk = (/ 2.293532201052922e-02_dp,   6.309209262997855e-02_dp, &
    !             1.047900103222502e-01_dp,   1.406532597155259e-01_dp, &
    !             1.690047266392679e-01_dp,   1.903505780647854e-01_dp, &
    !             2.044329400752989e-01_dp,   2.094821410847278e-01_dp    /)
    !    wg  = (/ 1.294849661688697e-01_dp,   2.797053914892767e-01_dp, &
    !             3.818300505051189e-01_dp,   4.179591836734694e-01_dp    /)

    !    centr  = 5.0e-01_dp*(a+b)
    !    hlgth  = 5.0e-01_dp*(b-a)
    !    dhlgth = abs(hlgth)

    !    ! Compute the 15-point Kronrod approximation to the integral,
    !    ! and estimate the absolute error.
    !    fc = f(centr)
    !    resg = fc*wg(4)
    !    resk = fc*wgk(8)
    !    resabs = abs(resk)

    !    do j = 1, 3
    !       jtw = j*2
    !       absc = hlgth*xgk(jtw)
    !       fval1 = f(centr-absc)
    !       fval2 = f(centr+absc)
    !       fv1(jtw) = fval1
    !       fv2(jtw) = fval2
    !       fsum = fval1+fval2
    !       resg = resg+wg(j)*fsum
    !       resk = resk+wgk(jtw)*fsum
    !       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    !    enddo

    !    do j = 1, 4
    !       jtwm1 = j*2-1
    !       absc = hlgth*xgk(jtwm1)
    !       fval1 = f(centr-absc)
    !       fval2 = f(centr+absc)
    !       fv1(jtwm1) = fval1
    !       fv2(jtwm1) = fval2
    !       fsum = fval1+fval2
    !       resk = resk+wgk(jtwm1)*fsum
    !       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    !    enddo

    !    reskh = resk * 5.0e-01_dp
    !    resasc = wgk(8)*abs(fc-reskh)

    !    do j = 1, 7
    !       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    !    enddo

    !    result = resk*hlgth
    !    resabs = resabs*dhlgth
    !    resasc = resasc*dhlgth
    !    abserr = abs((resk-resg)*hlgth)

    !    if ( resasc /= 0.0_dp.and.abserr /= 0.0_dp ) &
    !       abserr = resasc*min ( 1.0_dp,(2.0e+02_dp*abserr/resasc)**1.5_dp)

    !    if ( resabs > tiny ( resabs ) / (5.0e+01_dp* epsilon ( resabs ) ) ) &
    !       abserr = max (( epsilon ( resabs ) *5.0e+01_dp)*resabs,abserr)

    !end subroutine qk15


    ! ==========================================================================
    ! qags routine from quadpack
    ! Adaptive integration is used for the GPD-independent integrals at NLO
    ! in tiktaalik

    function iqags(f,a,b,oepsabs,oepsrel) result(res)
        ! Easy interface to qags.
        ! Written by Ian Cloet sometime before 2017.
        real(dp), external             :: f
        real(dp), intent(in)           :: a,b
        real(dp), intent(in), optional :: oepsabs,oepsrel
        integer  :: number_evaluations,ifail
        real(dp) :: res,epsabs,epsrel,absolute_error
        epsabs = 1.0e-12_dp
        epsrel = 1.0e-06_dp
        if(present(oepsabs)) epsabs = oepsabs
        if(present(oepsrel)) epsrel = oepsrel
        call qags(f,a,b,epsabs,epsrel,res,absolute_error,number_evaluations,ifail)
    end function iqags

    recursive subroutine qags(f,a,b,epsabs,epsrel,result,abserr,oneval,oifail)
        ! QAGS estimates the integral of a function.
        !
        ! Discussion:
        !
        !   The routine calculates an approximation RESULT to a definite integral
        !     I = integral of F over (A,B),
        !   hopefully satisfying
        !     || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
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
        !   Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
        !
        !   Output, real RESULT, the estimated value of the integral.
        !
        !   Output, real ABSERR, an estimate of || I - RESULT ||.
        !
        !   Output, integer NEVAL, the number of times the integral was evaluated.
        !
        !   Output, integer ifAIL, error flag.
        !                  ifail = 0 normal and reliable termination of the
        !                            routine. it is assumed that the requested
        !                            accuracy has been achieved.
        !                  ifail > 0 abnormal termination of the routine
        !                            the estimates for integral and error are
        !                            less reliable. it is assumed that the
        !                            requested accuracy has not been achieved.
        !                        = 1 maximum number of subdivisions allowed
        !                            has been achieved. one can allow more sub-
        !                            divisions by increasing the data value of
        !                            limit in qags (and taking the according
        !                            dimension adjustments into account).
        !                            however, if this yields no improvement
        !                            it is advised to analyze the integrand
        !                            in order to determine the integration
        !                            difficulties. if the position of a
        !                            local difficulty can be determined (e.g.
        !                            singularity, discontinuity within the
        !                            interval) one will probably gain from
        !                            splitting up the interval at this point
        !                            and calling the integrator on the sub-
        !                            ranges. if possible, an appropriate
        !                            special-purpose integrator should be used,
        !                            which is designed for handling the type
        !                            of difficulty involved.
        !                        = 2 the occurrence of roundoff error is detec-
        !                            ted, which prevents the requested
        !                            tolerance from being achieved.
        !                            the error may be under-estimated.
        !                        = 3 extremely bad integrand behavior occurs
        !                            at some  points of the integration
        !                            interval.
        !                        = 4 the algorithm does not converge. roundoff
        !                            error is detected in the extrapolation
        !                            table. it is presumed that the requested
        !                            tolerance cannot be achieved, and that the
        !                            returned result is the best which can be
        !                            obtained.
        !                        = 5 the integral is probably divergent, or
        !                            slowly convergent. it must be noted that
        !                            divergence can occur with any other value
        !                            of ifail.
        !                        = 6 the input is invalid, because
        !                            epsabs < 0 and epsrel < 0,
        !                            result, abserr and neval are set to zero.
        !
        ! Local Parameters:
        !
        !          alist     - list of left end points of all subintervals
        !                      considered up to now
        !          blist     - list of right end points of all subintervals
        !                      considered up to now
        !          rlist(i)  - approximation to the integral over
        !                      (alist(i),blist(i))
        !          rlist2    - array of dimension at least limexp+2 containing
        !                      the part of the epsilon table which is still
        !                      needed for further computations
        !          elist(i)  - error estimate applying to rlist(i)
        !          maxerr    - pointer to the interval with largest error
        !                      estimate
        !          errmax    - elist(maxerr)
        !          erlast    - error on the interval currently subdivided
        !                      (before that subdivision has taken place)
        !          area      - sum of the integrals over the subintervals
        !          errsum    - sum of the errors over the subintervals
        !          errbnd    - requested accuracy max(epsabs,epsrel*
        !                      abs(result))
        !          *****1    - variable for the left interval
        !          *****2    - variable for the right interval
        !          last      - index for subdivision
        !          nres      - number of calls to the extrapolation routine
        !          numrl2    - number of elements currently in rlist2. if an
        !                      appropriate approximation to the compounded
        !                      integral has been obtained it is put in
        !                      rlist2(numrl2) after numrl2 has been increased
        !                      by one.
        !          small     - length of the smallest interval considered
        !                      up to now, multiplied by 1.5
        !          erlarg    - sum of the errors over the intervals larger
        !                      than the smallest interval considered up to now
        !          extrap    - logical variable denoting that the routine is
        !                      attempting to perform extrapolation i.e. before
        !                      subdividing the smallest interval we try to
        !                      decrease the value of erlarg.
        !          noext     - logical variable denoting that extrapolation
        !                      is no longer allowed (true value)

        real(dp), external              :: f
        real(dp), intent(in)            :: a,b,epsabs,epsrel
        real(dp), intent(out)           :: result,abserr
        integer,  intent(out), optional :: oneval,oifail


        integer, parameter :: limit = 500
        logical  :: extrap,noext
        integer  :: id,ierro,iord(limit),iroff1,iroff2,iroff3,jupbnd,k,ksgn,ktmin, &
                  & last,maxerr,nres,nrmax,numrl2,neval,ifail
        real(dp) :: abseps,alist(limit),area,area1,area12,area2,a1,a2,blist(limit),b1,b2, &
                  & correc,defabs,defab1,defab2,dres,elist(limit),erlarg,erlast,errbnd,   &
                  & errmax,error1,error2,erro12,errsum,ertest,resabs,reseps,res3la(3),    &
                  & rlist(limit),rlist2(52),small

        ! The dimension of rlist2 is determined by the value of
        ! limexp in QEXTR (rlist2 should be of dimension (limexp+2) at least).

        ! Test on validity of parameters.
        ifail    = 0
        neval    = 0
        last     = 0
        result   = 0.0_dp
        abserr   = 0.0_dp
        alist(1) = a
        blist(1) = b
        rlist(1) = 0.0_dp
        elist(1) = 0.0_dp

        if ( epsabs < 0.0_dp .and. epsrel < 0.0_dp ) then
          ifail = 6
          return
        endif

        ! First approximation to the integral.
        ierro = 0
        call qk21(f,a,b,result,abserr,defabs,resabs)

        ! Test on accuracy.
        dres = abs ( result )
        errbnd = max ( epsabs, epsrel * dres )
        last = 1
        rlist(1) = result
        elist(1) = abserr
        iord(1) = 1

        if ( abserr <= 1.0e+02_dp * epsilon ( defabs ) * defabs .and.  abserr > errbnd ) ifail = 2
        if ( limit == 1 ) ifail = 1
        if ( ifail /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. abserr == 0.0_dp ) goto 140

        ! Initialization.
        rlist2(1) = result
        errmax = abserr
        maxerr = 1
        area   = result
        errsum = abserr
        abserr = huge ( abserr )
        nrmax  = 1
        nres   = 0
        numrl2 = 2
        ktmin  = 0
        extrap = .false.
        noext  = .false.
        iroff1 = 0
        iroff2 = 0
        iroff3 = 0
        correc = 0 ! Added by AF, Dec 11, 2024
        erlarg = 0 ! Added by AF, Dec 11, 2024
        small  = 0 ! Added by AF, Dec 11, 2024
        ertest = 0 ! Added by AF, Dec 11, 2024

        if ( dres >= (1.0_dp-5.0e+01_dp* epsilon ( defabs ) )*defabs ) then
          ksgn = 1
        else
          ksgn = -1
        endif

        do last = 2, limit ! AF: the corresponding enddo is on the line I marked "boop"

          ! Bisect the subinterval with the nrmax-th largest error estimate.
          a1 = alist(maxerr)
          b1 = 5.0e-01_dp*(alist(maxerr)+blist(maxerr))
          a2 = b1
          b2 = blist(maxerr)
          erlast = errmax
          call qk21(f,a1,b1,area1,error1,resabs,defab1)
          call qk21(f,a2,b2,area2,error2,resabs,defab2)

          ! Improve previous approximations to integral and error and test for accuracy.
          area12 = area1+area2
          erro12 = error1+error2
          errsum = errsum+erro12-errmax
          area = area+area12-rlist(maxerr)

          if ( defab1 == error1 .or. defab2 == error2 ) goto 15

          if ( abs ( rlist(maxerr) - area12) > 1.0e-05 * abs(area12) .or. erro12 < 9.9e-01_dp * errmax ) goto 10

          if ( extrap ) then
            iroff2 = iroff2+1
          else
            iroff1 = iroff1+1
          endif

10        continue
          if ( last > 10.and.erro12 > errmax ) iroff3 = iroff3+1

15        continue
          rlist(maxerr) = area1
          rlist(last) = area2
          errbnd = max ( epsabs,epsrel*abs(area))

          ! Test for roundoff error and eventually set error flag.
          if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) ifail = 2
          if ( iroff2 >= 5 ) ierro = 3

          ! Set error flag in the case that the number of subintervals equals limit.
          if ( last == limit ) ifail = 1

          ! Set error flag in the case of bad integrand behavior at a point of the integration range.
          if ( max(abs(a1),abs(b2)) <= (1.0_dp+1.0e+03_dp*epsilon(a1))*(abs(a2)+1.0e+03_dp*tiny(a2)) ) ifail = 4

          ! Append the newly-created intervals to the list.
          if ( error2 <= error1 ) then
            alist(last) = a2
            blist(maxerr) = b1
            blist(last) = b2
            elist(maxerr) = error1
            elist(last) = error2
          else
            alist(maxerr) = a2
            alist(last) = a1
            blist(last) = b1
            rlist(maxerr) = area2
            rlist(last) = area1
            elist(maxerr) = error2
            elist(last) = error1
          endif

          ! Call QSORT to maintain the descending ordering
          ! in the list of error estimates and select the subinterval
          ! with nrmax-th largest error estimate (to be bisected next).
          call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)
          if ( errsum <= errbnd ) goto 115 ! AF: escapes loop, computes integral
          if ( ifail /= 0 ) exit
          if ( last == 2 ) goto 80
          if ( noext ) goto 90
          erlarg = erlarg-erlast
          if ( abs(b1-a1) > small ) erlarg = erlarg+erro12
          if ( extrap ) goto 40

          ! Test whether the interval to be bisected next is the smallest interval.
          if ( abs(blist(maxerr)-alist(maxerr)) > small ) goto 90
          extrap = .true.
          nrmax = 2

40        continue
          ! The smallest interval has the largest error.
          ! Before bisecting decrease the sum of the errors over the
          ! larger intervals (erlarg) and perform extrapolation.
          if ( ierro /= 3 .and. erlarg > ertest ) then
            id = nrmax
            jupbnd = last
            if ( last > (2+limit/2) ) then
              jupbnd = limit+3-last
            endif
            do k = id, jupbnd
              maxerr = iord(nrmax)
              errmax = elist(maxerr)
              if ( abs(blist(maxerr)-alist(maxerr)) > small ) goto 90
              nrmax = nrmax+1
            enddo
          endif

          ! Perform extrapolation.
          ! AF: why is this line numbered? There are no goto statements to this.
60        continue
          numrl2 = numrl2+1
          rlist2(numrl2) = area
          call qextr(numrl2,rlist2,reseps,abseps,res3la,nres)
          ktmin = ktmin+1
          if ( ktmin > 5 .and. abserr < 1.0e-03_dp * errsum ) ifail = 5
          if ( abseps < abserr ) then
            ktmin  = 0
            abserr = abseps
            result = reseps
            correc = erlarg
            ertest = max ( epsabs,epsrel*abs(reseps))
            if ( abserr <= ertest ) exit
          endif

          ! Prepare bisection of the smallest interval.
          if ( numrl2 == 1 ) noext = .true.
          if ( ifail == 5 ) exit
          maxerr = iord(1)
          errmax = elist(maxerr)
          nrmax  = 1
          extrap = .false.
          small  = small*5.0e-01_dp
          erlarg = errsum
          goto 90

80        continue
          small  = abs(b-a)*3.75e-01_dp
          erlarg = errsum
          ertest = errbnd
          rlist2(2) = area

90        continue
        enddo ! boop ... end of the "do last = 2, limit" loop

        ! Set final result and error estimate.
        if ( abserr == huge ( abserr ) ) goto 115
        if ( ifail+ierro == 0 ) goto 110
        if ( ierro == 3 ) abserr = abserr+correc
        if ( ifail == 0 ) ifail = 3
        if ( result /= 0.0_dp.and.area /= 0.0_dp ) goto 105
        if ( abserr > errsum ) goto 115
        if ( area == 0.0_dp ) goto 130
        goto 110

105     continue
        if ( abserr/abs(result) > errsum/abs(area) ) goto 115

        ! Test on divergence.
110     continue
        if ( ksgn == (-1) .and. max(abs(result),abs(area)) <= defabs*1.0e-02_dp ) goto 130
        if ( 1.0e-02_dp > (result/area) .or. (result/area) > 1.0e+02_dp .or. errsum > abs(area) ) ifail = 6
        goto 130

        ! Compute global integral sum.
115     continue
        result = sum ( rlist(1:last) )
        abserr = errsum

130     continue
        if ( ifail > 2 ) ifail = ifail-1

140     continue
        neval = 42*last-21
        if (present(oneval)) oneval = neval
        if (present(oifail)) oifail = ifail

    end subroutine qags

    ! ==========================================================================
    ! Auxiliary routines from quadpack needed by qags

    subroutine qsort(limit,last,maxerr,ermax,elist,iord,nrmax)
        ! QSORT maintains the order of a list of local error estimates.
        !
        ! Discussion:
        !
        !   This routine maintains the descending ordering in the list of the
        !   local error estimates resulting from the interval subdivision process.
        !   At each call two error estimates are inserted using the sequential
        !   search top-down for the largest error estimate and bottom-up for the
        !   smallest error estimate.
        !
        ! Reference:
        !
        !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
        !   QUADPACK, a Subroutine Package for Automatic Integration,
        !   Springer Verlag, 1983
        !
        ! Parameters:
        !
        !   Input, integer LIMIT, the maximum number of error estimates the list can
        !   contain.
        !
        !   Input, integer LAST, the current number of error estimates.
        !
        !   Input/output, integer MAXERR, the index in the list of the NRMAX-th
        !   largest error.
        !
        !   Output, real ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
        !
        !   Input, real ELIST(LIMIT), contains the error estimates.
        !
        !   Input/output, integer IORD(LAST).  The first K elements contain
        !   pointers to the error estimates such that ELIST(IORD(1)) through
        !   ELIST(IORD(K)) form a decreasing sequence, with
        !     K = LAST
        !   if
        !     LAST <= (LIMIT/2+2),
        !   and otherwise
        !     K = LIMIT+1-LAST.
        !
        !   Input/output, integer NRMAX.

        integer,  intent(in)    :: limit,last
        real(dp), intent(in)    :: elist(limit)
        integer,  intent(inout) :: maxerr,nrmax,iord(last)
        real(dp), intent(out)   :: ermax

        ! Declare local variables
        integer  :: i,ibeg,isucc,j,jbnd,jupbn,k
        real(dp) :: errmax,errmin

        ! Check whether the list contains more than two error estimates.
        if ( last <= 2 ) then
          iord(1) = 1
          iord(2) = 2
          goto 90
        endif

        ! This part of the routine is only executed if, due to a
        ! difficult integrand, subdivision increased the error
        ! estimate. in the normal case the insert procedure should
        ! start after the nrmax-th largest error estimate.
        errmax = elist(maxerr)

        do i = 1, nrmax-1
          isucc = iord(nrmax-1)
          if ( errmax <= elist(isucc) ) exit
          iord(nrmax) = isucc
          nrmax = nrmax-1
        enddo

        ! Compute the number of elements in the list to be maintained
        ! in descending order.  This number depends on the number of
        ! subdivisions still allowed.
        jupbn = last
        if ( last > (limit/2+2) ) jupbn = limit+3-last
        errmin = elist(last)

        ! Insert errmax by traversing the list top-down, starting
        ! comparison from the element elist(iord(nrmax+1)).
        jbnd = jupbn-1
        ibeg = nrmax+1

        do i = ibeg, jbnd
          isucc = iord(i)
          if ( errmax >= elist(isucc) ) goto 60
          iord(i-1) = isucc
        enddo

        iord(jbnd) = maxerr
        iord(jupbn) = last
        goto 90

        ! Insert errmin by traversing the list bottom-up.
60      continue
        iord(i-1) = maxerr
        k = jbnd
        do j = i, jbnd
          isucc = iord(k)
          if ( errmin < elist(isucc) ) goto 80
          iord(k+1) = isucc
          k = k-1
        enddo
        iord(i) = last
        goto 90

80      continue
        iord(k+1) = last

        ! Set maxerr and ermax.
90      continue
        maxerr = iord(nrmax)
        ermax  = elist(maxerr)
    end subroutine qsort


    subroutine qextr(n,epstab,result,abserr,res3la,nres )
        ! QEXTR carries out the Epsilon extrapolation algorithm.
        !
        ! Discussion:
        !
        !   The routine determines the limit of a given sequence of approximations,
        !   by means of the epsilon algorithm of P. Wynn.  An estimate of the
        !   absolute error is also given.  The condensed epsilon table is computed.
        !   Only those elements needed for the computation of the next diagonal
        !   are preserved.
        !
        ! Reference:
        !
        !   R Piessens, E de Doncker-Kapenger, C W Ueberhuber, D K Kahaner,
        !   QUADPACK, a Subroutine Package for Automatic Integration,
        !   Springer Verlag, 1983
        !
        ! Parameters:
        !
        !   Input, integer N, indicates the entry of EPSTAB which contains
        !   the new element in the first column of the epsilon table.
        !
        !   Input/output, real EPSTAB(52), the two lower diagonals of the triangular
        !   epsilon table.  The elements are numbered starting at the right-hand
        !   corner of the triangle.
        !
        !   Output, real RESULT, the estimated value of the integral.
        !
        !   Output, real ABSERR, estimate of the absolute error computed from
        !   RESULT and the 3 previous results.
        !
        !   ?, real RES3LA(3), the last 3 results.
        !
        !   Input/output, integer NRES, the number of calls to the routine.  This
        !   should be zero on the first call, and is automatically updated
        !   before return.
        !
        ! Local Parameters:
        !
        !          e0     - the 4 elements on which the
        !          e1       computation of a new element in
        !          e2       the epsilon table is based
        !          e3                 e0
        !                       e3    e1    new
        !                             e2
        !          newelm - number of elements to be computed in the new
        !                   diagonal
        !          error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
        !          result - the element in the new diagonal with least value
        !                   of error
        !          limexp is the maximum number of elements the epsilon table
        !          can contain. if this number is reached, the upper diagonal
        !          of the epsilon table is deleted.

        integer,  intent(inout) :: n
        real(dp), intent(inout) :: epstab(52)
        integer,  intent(inout)   :: nres ! AF, Dec 11, 2024: was intent(out), fixed to intent(inout)
        real(dp), intent(out)   :: result,abserr,res3la(3)

        ! Declare local variables
        integer  :: i,ib,ib2,ie,indx,k1,k2,k3,limexp,newelm,num
        real(dp) :: delta1,delta2,delta3,epsinf,error,err1,err2,err3,e0,e1,e1abs,e2, &
                  & e3,res,ss,tol1,tol2,tol3

        nres   = nres+1
        abserr = huge ( abserr )
        result = epstab(n)

        if ( n < 3 ) goto 100
        limexp = 50
        epstab(n+2) = epstab(n)
        newelm = (n-1)/2
        epstab(n) = huge ( epstab(n) )
        num = n
        k1 = n

        do i = 1, newelm
          k2 = k1-1
          k3 = k1-2
          res = epstab(k1+2)
          e0 = epstab(k3)
          e1 = epstab(k2)
          e2 = res
          e1abs = abs(e1)
          delta2 = e2-e1
          err2 = abs(delta2)
          tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
          delta3 = e1-e0
          err3 = abs(delta3)
          tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )

          ! If e0, e1 and e2 are equal to within machine accuracy, convergence
          ! is assumed.
          if ( err2 <= tol2 .and. err3 <= tol3 ) then
            result = res
            abserr = err2+err3
            goto 100 ! AF: breaks out of loop
          endif

          e3 = epstab(k1)
          epstab(k1) = e1
          delta1 = e1-e3
          err1 = abs(delta1)
          tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )

          ! If two elements are very close to each other, omit a part
          ! of the table by adjusting the value of N.
          if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) goto 20
          ss = 1.0_dp/delta1+1.0_dp/delta2-1.0_dp/delta3
          epsinf = abs ( ss*e1 )

          ! Test to detect irregular behavior in the table, and
          ! eventually omit a part of the table adjusting the value of N.
          if ( epsinf > 1.0e-04_dp ) goto 30

20        continue
          n = i+i-1
          exit

          ! Compute a new element and eventually adjust the value of RESULT.
30        continue
          res = e1+1.0_dp/ss
          epstab(k1) = res
          k1 = k1-2
          error = err2+abs(res-e2)+err3
          if ( error <= abserr ) then
            abserr = error
            result = res
          endif

        enddo ! ends "do i = 1, newelm"

        ! Shift the table.
        if ( n == limexp ) then
          n = 2*(limexp/2)-1
        endif

        if ( (num/2)*2 == num ) then
          ib = 2
        else
          ib = 1
        endif

        ie = newelm+1

        do i = 1, ie
          ib2 = ib+2
          epstab(ib) = epstab(ib2)
          ib = ib2
        enddo

        if ( num /= n ) then
          indx = num-n+1
          do i = 1, n
            epstab(i)= epstab(indx)
            indx = indx+1
          enddo
        endif

        if ( nres < 4 ) then
          res3la(nres) = result
          abserr = huge ( abserr )
        else
          abserr = abs(result-res3la(3))+abs(result-res3la(2)) + abs(result-res3la(1))
          res3la(1) = res3la(2)
          res3la(2) = res3la(3)
          res3la(3) = result
        endif

100     continue
        abserr = max(abserr, 0.5_dp*epsilon(result)*abs(result))

    end subroutine qextr


end module integration

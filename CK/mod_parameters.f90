module mod_parameters


    ! 0: optimized, _ini: initial value to be updated, ?: flag for initial value
    !******* metaparameters *******!
    integer, parameter :: maxiter = 2000    ! the maximum number of value function iterations in innerloop
    integer, parameter :: diagnum = 500     ! display result in every diagnum interval (innerloop)
    integer, parameter :: diagnumout = 500  ! display result in every diagnum interval (outerloop)
    character(len=*), parameter :: jsonfilename = './results_CKnotr_KS.json'
    logical, parameter :: jsonoutput = .true.

    real(8), parameter :: tolout = 1d-4      ! tolerance level for market clearing in steady state
    real(8), parameter :: tol = 1d-7        ! tolerance level for innerloop and (tol^2) for optimization
    real(8), parameter :: toldist = 1d-10 !1d-7 !1d-5 !1d-7    ! tolerance level for distribution (for stationary equilibrium)
    ! real(8), parameter :: tolout = 1d-4     ! tolerance level for outer loop
    real(8), parameter :: tolmkt = 1d-5 !1d-6     ! tolerance level for market clearing in outerloop
    real(8), parameter :: dampss = 0.95d0  ! dampening weight for market clearing in stationary equilibrium (for old) note: values close to 0.5 (e.g., 0.7) lead to divergence
    real(8), parameter :: dampout = 0.7d0  !0.7d0  ! dampening weight for outer loop (for old)

    integer, parameter :: transmat = 3      ! =0 iteration method without transition matrix
                                            ! =1 iteration method with transition matrix
                                            ! =2 eigenvalue decomposition by lapack (with spblas=0)
                                            ! =3 eigenvalue decomposition by splapack (with spblas=1); need arpack library
    integer, parameter :: spblas = 1        ! =1 use sparse matrix (with transmat>=1)

    integer, parameter :: linflag = 2       ! =0 use 2 dim spline over (a,m),
                                            ! =1 use linear interpolation over (a,m),
                                            ! =2 use 1 dim spline over a and linear interpolation over m
    integer, parameter :: howflag = 1       ! =1 use Howard's improvement algorithm in innerloop
    ! for convergence criteria (both in the steady state and the inner loop)
    ! NOTE: if diff < 1, log(diff) (percantage difference) has more severe criteria as the slope of log(diff) gets steeper
    integer, parameter :: logdiff = 0       ! =1 percetage difference

    integer, parameter :: conzflag = 1      ! =3 use continuous z (linear interpolation) and ln(z) & ln(z)*ln(K) as regressors in forecasting rules
                                            ! =2 use continuous z (linear interpolation) and ln(z)-as a regressor in forecasting rules
                                            ! =1 use z-conditional forecasting rules
                                            ! =0 use ln(z) as a regressor in forecasting rules
    integer, parameter :: calcma = 1        ! =1 use minasset, the threshold below which hshlds choose to work
    integer, parameter :: outermkt = 0      ! =1 market clearing in outerloop
    integer, parameter :: maxcountmkt = 30  ! the maximum number of bisection iterations
    integer, parameter :: bsctr = 0         ! =1 bisection for market clearing over r, =0 over w
    integer, parameter :: bsfixr = 0        ! =1 fix r (w if bsctr=1) while updating w in market clearing
    real(8), parameter :: bspct = 0.01d0    ! initial range for bisection (the previous error is used after time 2)

    ! for the initial guess of the value function in the inner loop
    ! integer, parameter :: ssini0 = 0        ! =1 reading the previous steady state results
    integer, parameter :: vmatini = 1       ! =1 using the one in the previous iteration (used in innerloop)
    integer, parameter :: vmatini0 = 0      ! =1 using the one saved in Vmat.txt for the first iteration
    ! ! for the initial guess of the forecasting rules in the outer loop
    ! ! NOTE: 030318 The variable name below is changed
    integer, parameter :: fcstini0 = 0      ! =1 reading the forecasting rules in the files below
    character(40) :: ckappakp = "./kappakp.txt"
    character(40) :: ckappaw  = "./kappaw.txt"
    character(40) :: ckappar  = "./kappar.txt"
    character(40) :: ckappal  = "./kappal.txt"
    ! forecasting rules chosen to be estimated
    integer, parameter :: fcsteqn = 1
    ! =1: kappakp and kappaw (kappar is obtained by the firm's FOC)
    ! =2: kappakp and kappar (kappaw is obtained by the firm's FOC)
    ! =3: kappakp, kappaw, and kappar

    ! for degree adjustment
    integer, parameter :: congrid = 0       ! grid option       ! =1 custermized
    real(8), parameter :: deggrid = 4.0d0   ! concavity degree adjustment
    ! Approximation of individual productivity
    integer, parameter :: nx = 17 !17
    integer, parameter :: tauflag = 0 ! =1 Tauchen approximation, =0 Rouwenhorst instead
    real(8), parameter :: mx = 3.0d0 !3.0d0
    ! ! Approximation of aggregate productivity
    integer, parameter :: nz = 7
    integer, parameter :: tauflagagg = 0
    real(8), parameter :: mz = 3.0d0
    ! integer, parameter :: nghz = 5

    !******* the number of grid points and the length of simulation *******!
    integer, parameter :: nb = 1
    integer, parameter :: ne = nb*nx
    integer, parameter :: na = 50 !50 !100 !200             ! for value functions and policy functions
    integer, parameter :: nk = 500 !500 !1000 !2000            ! for distributions
    integer, parameter :: nm = 7
    real(8), parameter :: mpct = 0.15d0        ! percent range of aggregate capital
    integer, parameter :: ra = na-2
    integer, parameter :: rm = nm-2

    ! for simulation
    integer, parameter :: simT = 3000
    integer, parameter :: drop = 500
    integer, parameter :: simTT = simT+drop
    ! ! for impulse reponse analysis
    ! real(8), parameter :: shocksize = -0.02d0 ! 1% of TFP shock
    ! integer, parameter :: irT = 200
    ! integer, parameter :: irdrop = 100
    ! integer, parameter :: irTT = irT+irdrop

    !*******  structural parameters *******!
    real(8), parameter :: ALPHA = 0.36d0        ! Capital share
    real(8), parameter :: DELTA = 0.025d0       ! Depreciation rate
    real(8), parameter :: HBAR = 1.0d0/3.0d0    ! hours if working
    real(8), parameter :: RHOZ = 0.95d0         ! Persistence of aggregate shocks
    real(8), parameter :: SDINOVZ = 0.007d0     ! SD of innovations to aggregate shocks
    ! Tax
    real(8), parameter :: taul = 0.279d0

    ! Calibrated using SMM
    ! Model A: Baseline with productivity dependent specification (OK)
    ! NOTE: With the parameters below and the transfer function, the steady state code doesn't converge
    ! from August 2018 version
    ! real(8), parameter :: B0 = 0.8760d0
    ! real(8), parameter :: BETA = 0.9863d0
    ! real(8), parameter :: RHOX = 0.9816d0         ! Persistence of idiosyncratic shocks
    ! real(8), parameter :: SDINOVX = 0.0931d0     ! SD of innovations to idiosyncratic shocks
    ! real(8), parameter :: T0 = 0.0705d0
    ! real(8), parameter :: Trscale = 0.498d0
    ! real(8), parameter :: Transprog = 4.08d0
    ! real(8), parameter :: phi = 0.0d0

    ! Model B: CK version (OK)
    ! HA-N: 20181125 New targets (age 23-70)
    real(8), parameter :: B0 = 0.880d0
    real(8), parameter :: BETA = 0.9848d0
    real(8), parameter :: RHOX = 0.9757d0        ! Persistence of idiosyncratic shocks
    real(8), parameter :: SDINOVX = 0.132d0     ! SD of innovations to idiosyncratic shocks
    ! August
    ! real(8), parameter :: B0 = 1.092d0
    ! real(8), parameter :: BETA = 0.9841d0
    ! real(8), parameter :: RHOX = 0.9726d0        ! Persistence of idiosyncratic shocks
    ! real(8), parameter :: SDINOVX = 0.1310d0     ! SD of innovations to idiosyncratic shocks
    ! May
    ! real(8), parameter :: B0 = 1.041d0
    ! real(8), parameter :: BETA = 0.9829d0
    ! real(8), parameter :: RHOX = 0.9719d0        ! Persistence of idiosyncratic shocks
    ! real(8), parameter :: SDINOVX = 0.156d0     ! SD of innovations to idiosyncratic shocks
    real(8), parameter :: Trscale = 0.0d0
    real(8), parameter :: Transprog = 0.0d0
    real(8), parameter :: T0 = 0.0d0
    real(8), parameter :: phi = 0.0d0


    ! Note: not used  ******************
    ! Estimates for the gs type progressive tax function in Guner and Venture RED (2014)
    ! a0[y-(y^-a1 + a2)^(-(1/a1))]
    real(8), parameter :: a0 = 0.276d0
    real(8), parameter :: a1 = 0.927d0
    real(8), parameter :: a2 = 4.15d0
    real(8), parameter :: wgtT = 0.0d0


end module mod_parameters

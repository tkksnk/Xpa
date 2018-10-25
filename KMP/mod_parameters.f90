 module mod_parameters

    integer, parameter :: maxiter = 2000
    integer, parameter :: diagnum = 200
    integer, parameter :: diagnumout = 500
    character(len=*), parameter :: jsonfilename = './results_calibKS_Reiter.json'
    ! character(len=*), parameter :: jsonfilename = './results_calibKMPbenchmark_Reiter.json'
    logical, parameter :: jsonoutput = .true.

    logical, parameter :: nraflag = .true.
    ! NOTE: 180801: linflag=1 and spliflag=1 doesn't work???
    logical, parameter :: linflag = .false.
    logical, parameter :: spliflag = .true. !.false.
    ! logical, parameter :: bisectm = .true.
    logical, parameter :: bisectm = .false.
    integer, parameter :: transmat = 1
    logical, parameter :: adjbias = .true.
    ! logical, parameter :: adjbias = .false.
    logical, parameter :: naiveflag = .false.
    integer, parameter :: congrid = 0
    real(8), parameter :: deggrid = 4.0d0
    real(8), parameter :: damp = 0.35d0 !0.35d0 ! for new: better convergence than 0.4
    real(8), parameter :: dampss = 0.005d0 !0.05d0 !0.005d0 ! for iterative method
    integer, parameter :: simT = 2000
    integer, parameter :: drop = 500
    integer, parameter :: simTT = simT+drop
    integer, parameter :: irfT = 100
    integer, parameter :: irfdrop = 2000
    integer, parameter :: irfTT = irfT+irfdrop

    integer, parameter :: nk = 25 !51 !101 !51 !25
    integer, parameter :: nm = 5
    integer, parameter :: nz = 2
    integer, parameter :: ne = 2
    integer, parameter :: ny = 1
    integer, parameter :: nd = 1
    ! real(8), parameter :: kmax = 20.0d0
    real(8), parameter :: kmax = 250.0d0 ! KS
    ! real(8), parameter :: kmax = 400.0d0 ! KS
    ! real(8), parameter :: kmax = 700.0d0 ! hety
    ! real(8), parameter :: kmax = 2000.0d0 ! bench
    integer, parameter :: nx = ne*ny*nd
    integer, parameter :: nb = 51 !201 !2001 ! NOTE: with KS calibration, kmax=200 and nb=5000 yields very different distribution
    integer, parameter :: rk = nk-2
    integer, parameter :: rm = nm-2

! parameters
    ! real(8), parameter :: THETA = 1.0d0-1.0d0/160.0d0
    real(8), parameter :: THETA = 1.0d0
    ! real(8), parameter :: BETA  = 0.96d0
    ! real(8), parameter :: SIGMA  = 1.5d0
    ! real(8), parameter :: ALPHA  = 0.36d0
    ! real(8), parameter :: DELTA  = 0.08d0
    real(8), parameter :: SIGMA  = 1.0d0
    real(8), parameter :: ALPHA  = 0.36d0
    real(8), parameter :: DELTA  = 0.025d0

    ! transition matrix for z
    ! ! KS 1998
    ! real(8), parameter :: zg = 1.01d0
    ! real(8), parameter :: zb = 0.99d0
    ! real(8), parameter :: ug = 0.04d0        ! unemployment rate at good time
    ! real(8), parameter :: ub = 0.10d0        ! unemployment rate at bad time
    ! real(8), parameter :: FractionZb = 0.5d0 ! fraction of time in severe recession
    ! real(8), parameter :: DurationZg = 8.0d0 ! average duration of good time
    ! real(8), parameter :: DurationZb = 8.0d0 ! average duration of bad time
    ! real(8), parameter :: pgg = 1.0d0-1.0d0/DurationZg
    ! real(8), parameter :: pbb = 1.0d0-1.0d0/DurationZb
    ! ! NOTE: hours worked in steady state h*(1-uss) is normalized to one
    ! real(8), parameter :: uss = (1.0d0-FractionZb)*ug + FractionZb*ub
    ! real(8), parameter :: mu = 0.15d0        ! unemployment insurance
    ! real(8), parameter :: h = 1.0d0/(1.0d0-uss)
    ! ! only those employed are taxed
    ! real(8), parameter :: taug = ug*mu/(h*(1.0d0-ug))
    ! real(8), parameter :: taub = ub*mu/(h*(1.0d0-ub))
    ! KMP 2015/17
    real(8), parameter :: zg = 1.0064d0
    real(8), parameter :: zb = 0.9676d0
    real(8), parameter :: ug = 1.0d0-0.946655772148944d0 !0.053344227851056d0      ! unemployment rate at good time
    real(8), parameter :: ub = 1.0d0-0.916159380188157d0 !0.083840619811843d0      ! unemployment rate at bad time
    real(8), parameter :: FractionZb = 0.1648d0 ! fraction of time in severe recession
    real(8), parameter :: DurationZb = 22.0d0 ! average duration of severe recession
    real(8), parameter :: pbb = 1.0d0-1.0d0/DurationZb
    real(8), parameter :: pgg = (1.0d0-FractionZb*(2.0d0-pbb))/(1.0d0-FractionZb)
    ! NOTE: hours worked in steady state h*(1-uss) is normalized to one
    real(8), parameter :: uss = (1.0d0-FractionZb)*ug + FractionZb*ub
    ! real(8), parameter :: mu = 0.5d0          ! unemployment insurance
    real(8), parameter :: mu = 0.01d0         ! unemployment insurance
    real(8), parameter :: h = 1.0d0/(1.0d0-ub) ! hours worked per worker
    real(8), parameter :: taug = ug*mu/(h*(1-ug)+ug*mu)
    real(8), parameter :: taub = ub*mu/(h*(1-ub)+ub*mu)

    ! transition matrix for e
    ! ! Aiyagari
    ! real(8), parameter :: mu = 0.05d0         ! unemployment insurance
    ! real(8), parameter :: pee = 0.925d0
    ! real(8), parameter :: puu = 0.5d0
    ! ! KS 1998
    ! real(8), parameter :: DurationUg = 1.5d0 ! average duration of unemployment at good time
    ! real(8), parameter :: DurationUb = 2.5d0 ! average duration of unemployment at bad time
    ! real(8), parameter :: corr = 0.25d0      ! correlation
    ! KMP 2015/17
    real(8), parameter :: pgg00 = 0.1890d0
    real(8), parameter :: pgg10 = 1.0d0-0.9543d0
    real(8), parameter :: pgb00 = 0.3382d0
    real(8), parameter :: pgb10 = 1.0d0-0.9304d0
    real(8), parameter :: pbg00 = 0.2220d0
    real(8), parameter :: pbg10 = 1.0d0-0.9622d0
    real(8), parameter :: pbb00 = 0.3378d0
    real(8), parameter :: pbb10 = 1.0d0-0.9394d0

    ! KMP 2015
    real(8), parameter :: rhoy = 0.9457d0
    real(8), parameter :: sigy = 0.0359d0**.5d0
    real(8), parameter :: dbar = 0.98349d0
    real(8), parameter :: deps = 0.01004d0
    ! ! KMP 2017?
    ! real(8), parameter :: rhoy = 0.9695d0
    ! real(8), parameter :: sigy = 0.0384**.5d0
    ! real(8), parameter :: dbar = 0.9864d0
    ! real(8), parameter :: deps = 0.0053d0

! inner/outer loop: critout>critin
    real(8), parameter :: critout = 1d-4
    real(8), parameter :: critbp  = 1d-5
    real(8), parameter :: critin  = 1d-6
    real(8), parameter :: critg   = 1d-10 ! machine precision?
    real(8), parameter :: critn   = 1d-10
    real(8), parameter :: critmu  = 1d-10 !1d-5

end module mod_parameters

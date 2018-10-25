 module mod_parameters

    integer, parameter :: diagnum = 200 ! for ss and inner loop
    integer, parameter :: diagnumout = 500 ! for outer loop
    ! character(len=*), parameter :: jsonfilename = './results_extend_KS.json'
    character(len=*), parameter :: jsonfilename = './results_extend_KS.json'
    logical, parameter :: jsonoutput = .true.

    logical, parameter :: nraflag = .false. ! NOTE: in the inner loop, nra sometime doesn't work... in planner, nra doesn't work
    logical, parameter :: bisectw = .true.
    logical, parameter :: outermkt = .true.
    logical, parameter :: fcstini = .false.
    logical, parameter :: adjbias = .true.
    ! logical, parameter :: adjbias = .false.
    ! logical, parameter :: naiveflag = .true.
    logical, parameter :: naiveflag = .false.
    real(8), parameter :: damp = 1.0d0
    real(8), parameter :: dampss = 0.05d0 ! for interation ss
    integer, parameter :: simT = 2000
    integer, parameter :: drop = 500
    integer, parameter :: simTT = simT+drop
    integer, parameter :: irfT = 50
    integer, parameter :: irfdrop = 200
    integer, parameter :: irfTT = irfT+irfdrop

    integer, parameter :: nk = 101
    integer, parameter :: nm = 5
    integer, parameter :: nz = 5
    integer, parameter :: ne = 5 !1
    integer, parameter :: nb = 2001 !2001
    integer, parameter :: rk = nk-2
    integer, parameter :: rm = nm-2
    real(8), parameter :: mz = 1.9600d0 ! 95%
    real(8), parameter :: me = 2.5758d0 ! 99%

! parameters, taken from p.406 of Khan and Thomas (2008)
    real(8), parameter :: GAMY = 1.016d0
    ! real(8), parameter :: GAMY = 1.0d0 !1.016d0
    real(8), parameter :: BETA = 0.97692d0
    real(8), parameter :: DELTA = 0.069d0
    real(8), parameter :: NU = 0.640d0
    real(8), parameter :: RHO = 0.85904793659574d0
    real(8), parameter :: SIGMA = 0.013961d0
    ! real(8), parameter :: THETA = 0.2568d0
    real(8), parameter :: THETA = 0.25648d0
    ! real(8), parameter :: ETA = 2.406d0
    real(8), parameter :: ETA = 2.400d0
    real(8), parameter :: RHOE = 0.85904793659574d0
    real(8), parameter :: SIGE = 0.022d0
    ! real(8), parameter :: RHOE = 0.0d0 !0.85904793659574d0
    ! real(8), parameter :: SIGE = 1d-6 !0.022d0
    ! if (ne==1) then
        ! real(8), parameter :: B = 1d-6
        ! real(8), parameter :: XIBAR = 0.014d0
    ! else
        real(8), parameter :: B = 0.011d0
        real(8), parameter :: XIBAR = 0.00825d0 ! p. 426
    ! end if
    real(8), parameter :: shocksize = SIGMA

! inner/outer loop: critout>critin
    real(8), parameter :: critout = 1d-4
    real(8), parameter :: critin  = 1d-5
    real(8), parameter :: critbp  = 1d-6
    real(8), parameter :: critbn  = 1d-6
    real(8), parameter :: critg   = 1d-10
    real(8), parameter :: critn   = 1d-10
    real(8), parameter :: critmu  = 1d-10 !1d-5

end module mod_parameters

module mod_planner


use mod_functions
use mod_spline
implicit none


contains


subroutine solvepl(knotsk,invTk,Gz,Pz,gmat0,nmat0,vmat0)


    use mod_parameters
    real(8), intent(in) :: knotsk(:), invTk(:,:), Gz(:), Pz(:,:)
    real(8), intent(out) :: gmat0(:,:), nmat0(:,:), vmat0(:,:)
    integer iz, jz, ik, iter, s1
    real(8) nLmat(nk,nz), nHmat(nk,nz), nmat1(nk,nz), vmat1(nk,nz), gmat1(nk,nz), vcond(nk), cf(4,rk+1), &
    know, znow, yterm, cnow, nnow, f0, df0, d2f0
    real(8) diffg, diffv, diff


    ! initial guess for vmat0 and nmat0
    do iz = 1,nz

        znow = Gz(iz)

        do ik = 1,nk

            know = knotsk(ik)

            nmat0(ik,iz) = bisectn(0.0d0,1.0d0,0.0d0,znow,know)
            nLmat(ik,iz) = bisectn(0.0d0,1.0d0,knotsk(1),znow,know)
            nHmat(ik,iz) = bisectn(0.0d0,1.0d0,knotsk(nk),znow,know)
            yterm = znow*know**THETA*nmat0(ik,iz)**NU + (1.0d0-DELTA)*know
            vmat0(ik,iz) = log(yterm) + ETA*(1.0d0-nmat0(ik,iz))

        end do

    end do


    gmat0 = 0.0d0

    diff = 1d+4
    iter = 0
    s1   = 0

    do while (diff>critin)

        !$omp parallel do private(znow,vcond,cf)
        do iz = 1,nz

            znow = Gz(iz)

            vcond = 0.0d0
            do jz = 1,nz

                vcond = vcond + Pz(iz,jz)*vmat0(:,jz)

            end do

            cf = spfit(invTk,vcond,rk,knotsk)

            !$omp parallel do private(know,nnow,yterm,cnow,f0,df0,d2f0)
            do ik = 1,nk

                know  = knotsk(ik)

                if (s1==0) then

                ! NOTE: Newton-Rhapson is to be done
                ! if (nraflag) then
                !     call nra(nmat0(ik,iz),nLmat(ik,iz),nHmat(ik,iz),nnow,cf,knotsk,know,znow)
                ! else
                    call gss(nLmat(ik,iz),nHmat(ik,iz),nnow,cf,knotsk,know,znow)
                ! end if

                    yterm = znow*know**THETA*nnow**NU + (1.0d0-DELTA)*know
                    cnow = (NU/ETA)*znow*know**THETA*nnow**(NU-1.0d0)

                    gmat1(ik,iz) = (yterm-cnow)/GAMY
                    nmat1(ik,iz) = nnow

                end if

                call vfuncsp(nmat1(ik,iz),cf,knotsk,know,znow,f0,df0,d2f0)
                vmat1(ik,iz) = -f0

            end do ! loop for k
            !$omp end parallel do

        end do ! loop for z
        !$omp end parallel do

        diffg = maxval(maxval(abs(gmat1-gmat0),1),1)
        diffv = maxval(maxval(abs(vmat1-vmat0),1),1)
        diff  = diffv
        iter = iter + 1
        ! diagnosis
        if (mod(iter,100)==0) then

            write(*,"('  iteration ', I4, '  ||Tv-v|| = ', F10.5, '  ||Tg-g|| = ', F10.5)") iter, diffv, diffg

        end if

        if (diffg<1d-4 .and. s1==0) then
            s1 = 1
        elseif (s1>0 .and. s1<20) then
            s1 = s1+1
        elseif (s1>=20) then
            s1 = 0
        end if

        nmat0 = nmat1
        gmat0 = gmat1
        vmat0 = vmat1

    end do


end subroutine solvepl


subroutine forecastpl(nmat0,knotsk,knotsm,Gz,mpmat0,pmat0)


    use mod_parameters
    real(8), intent(in) :: nmat0(:,:), knotsk(:), knotsm(:), Gz(:)
    real(8), intent(out) :: mpmat0(:,:), pmat0(:,:)
    real(8) znow, know, wk, nnow, yterm, cnow, kp
    integer iz, im, ik


    do iz = 1,nz

        znow = Gz(iz)

        do im = 1,nm

            know = knotsm(im)

            ! linear interpolation
            ik = gridlookup2(know,knotsk)
            ! ik = gridlookup(nk,knotsk,know)
            wk = (knotsk(ik+1)-know)/(knotsk(ik+1)-knotsk(ik))
            nnow = wk*nmat0(ik,iz) + (1.0d0-wk)*nmat0(ik+1,iz)

            yterm = znow*know**THETA*nnow**NU + (1.0d0-DELTA)*know

            ! first-order condition of n
            cnow = (NU/ETA)*znow*know**THETA*nnow**(NU-1.0d0)
            kp = (yterm-cnow)/GAMY

            mpmat0(im,iz) = kp
            pmat0(im,iz)  = 1.0d0/cnow

        end do

    end do


end subroutine forecastpl


subroutine simulatepl(nmat0,knotsk,Gz,izvec,Kvec,Kpvec,Cvec,Nvec)


    use mod_parameters
    real(8), intent(in)  :: nmat0(:,:), knotsk(:), Gz(:)
    integer, intent(in)  :: izvec(:)
    real(8), intent(out) :: Kvec(:), Kpvec(:), Cvec(:), Nvec(:)
    integer ik, iz, tt
    real(8) ykSS, ckSS, ycSS, nSS, kSS
    real(8) know, znow, wk, nnow, yterm, cnow, kp


    ! steady state
    ykSS = (GAMY-BETA*(1.0d0-DELTA))/BETA/THETA
    ckSS = ykSS + (1.0d0-GAMY-DELTA)
    ycSS = ykSS/ckSS
    nSS = NU/ETA*ycSS
    kSS = (ykSS*nSS**(-NU))**(1.0d0/(THETA-1.0d0))

    know = kSS

    do tt = 1,simTT !tott

        iz   = izvec(tt)
        znow = Gz(iz)

        ! linear interpolation
        ik = gridlookup2(know,knotsk)
        wk = (knotsk(ik+1)-know)/(knotsk(ik+1)-knotsk(ik))
        nnow = wk*nmat0(ik,iz) + (1.0d0-wk)*nmat0(ik+1,iz)

        yterm = znow*know**THETA*nnow**NU + (1.0d0-DELTA)*know

        ! first-order condition of n
        cnow = (NU/ETA)*znow*know**THETA*nnow**(NU-1.0d0)
        kp = (yterm-cnow)/GAMY

        ! record aggregate variables
        Kvec(tt) = know
        Kpvec(tt) = kp
        Cvec(tt) = cnow
        Nvec(tt) = nnow

        know = kp

    end do


end subroutine simulatepl


function bisectn(nlow,nhigh,kp,znow,know) result(n0)


    use mod_parameters, only: critbn
    real(8), intent(in) :: nlow, nhigh, kp, znow, know
    real(8) nL, nH, n0, B0, diff
    integer iterbn


    nL = nlow
    nH = nhigh
    diff = 1d+4
    iterbn = 0

    do while (diff>critbn)

        n0 = (nL+nH)/2.0d0
        B0 = nfunc(n0,kp,znow,know)

        if (B0<0.0d0) then
            nL = n0
        else
            nH = n0
        end if

        diff = nH-nL
        iterbn = iterbn + 1
    !     s = sprintf('  bisection %4d,  pH-pL = %6.8f',iterbn,diff);
    !     disp(s);

    end do


end function bisectn


function nfunc(n0,kp,znow,know) result(f)


    use mod_parameters
    real(8), intent(in) :: n0, kp, znow, know
    real(8) yterm, c, f


    yterm = znow*know**THETA*n0**NU + (1.0d0-DELTA)*know
    c = yterm - GAMY*kp

    ! first-order condition of n
    f = -NU*znow*know**THETA*n0**(NU-1.0d0) + ETA*c


end function nfunc


subroutine nra(xinit,xmin,xmax,x,cf,knotsk,know,znow)


    use mod_parameters, only: rk, critn
	real(8), intent(in) :: xinit, xmin, xmax, cf(:,:), knotsk(:), know, znow
	real(8), intent(out) :: x
    real(8) FL, DFL, FH, DFH, f0, df0, d2f0, NewtonStep, x0, xlow, xhigh, x1, diff
    integer iter


    call vfuncsp(xmin,cf,knotsk,know,znow,FL,DFL,d2f0)
    call vfuncsp(xmax,cf,knotsk,know,znow,FH,DFH,d2f0)

    if (DFL<0) then

        x0 = xmin

    elseif (DFH>0) then

        x0 = xmax

    else ! DFL>=0 amd DFH<=0

        x0    = xinit
        xlow  = xmin
        xhigh = xmax

        diff = 1d+4
        iter = 0

        ! Newton-Rhapson
        do while (diff>critn)

            call vfuncsp(x0,cf,knotsk,know,znow,f0,df0,d2f0)
            NewtonStep = df0/d2f0
            x1 = x0 - NewtonStep

            if (df0>0) then
                xlow = x0
            else
                xhigh = x0
            end if

            ! backtrack
            if (x1>xhigh) then
                x1 = (xlow+xhigh)/2.0d0
            elseif (x1<xlow) then
                x1 = (xlow+xhigh)/2.0d0
            end if

            diff = abs(x1-x0)
            iter = iter + 1
            x0 = x1

        end do

        ! check bounds
        if (x0<xmin) then
            x0 = xmin
        elseif (x0>xmax) then
            x0 = xmax
        end if

    end if

    x = x0


end subroutine nra


subroutine gss(xmin,xmax,x,cf,knotsk,know,znow)


    use mod_parameters, only: rk, critg
	real(8), intent(in) :: xmin, xmax, cf(:,:), knotsk(:), know, znow
	real(8), intent(out) :: x
	real(8) rg, a, b, c, d, z, fc, fd, crit, diff, df0, d2f0
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    call vfuncsp(c,cf,knotsk,know,znow,fc,df0,d2f0)
	d = a + (1.0d0-rg)*(b-a)
    call vfuncsp(d,cf,knotsk,know,znow,fd,df0,d2f0)

    diff = 1d+4
	iter = 0

	do while (diff>critg)

 		if (fc>=fd) then

        	z = c + (1.0d0-rg)*(b-c)
        	a = c
        	c = d
        	fc = fd
        	d = z
            call vfuncsp(d,cf,knotsk,know,znow,fd,df0,d2f0)

	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            call vfuncsp(c,cf,knotsk,know,znow,fc,df0,d2f0)

    	end if

	    diff = d-c
    	iter = iter+1

        ! diagnosis
!        write(*,"('  iteration ', I4, '  New x = ', F5.10)") iter, diff
        ! print *, c

	end do

!    pause

	x = c


end subroutine gss


subroutine vfuncsp(n0,cf,knotsk,know,znow,f0,df0,d2f0)


    use mod_parameters
	real(8), intent(in) :: n0, cf(4,rk+1), knotsk(nk), know, znow
    real(8), intent(out) :: f0, df0, d2f0
    real(8) yterm, cnow, kp, EV, EDV, ED2V, dc, d2c, dk, d2k

    ! first-order condition of n
    cnow = (NU/ETA)*znow*know**THETA*n0**(NU-1.0d0)
    dc   = (NU-1.0d0)*(NU/ETA)*znow*know**THETA*n0**(NU-2.0d0)
    d2c  = (NU-1.0d0)*(NU-2.0d0)*(NU/ETA)*znow*know**THETA*n0**(NU-3.0d0)

    yterm = znow*know**THETA*n0**NU + (1.0d0-DELTA)*know
    kp  = (yterm-cnow)/GAMY
    dk  = (NU*znow*know**THETA*n0**(NU-1.0d0)-dc)/GAMY
    d2k = (NU*(NU-1.0d0)*znow*know**THETA*n0**(NU-2.0d0)-d2c)/GAMY
    call speva(cf,kp,rk,knotsk,EV,EDV,ED2V)

    f0   = log(cnow) + ETA*(1.0d0-n0) + BETA*EV
    f0   = -f0
    df0  = 1.0d0/cnow*dc - ETA + BETA*EDV*dk
    d2f0 = -(1.0d0/cnow**2)*dc + 1.0d0/cnow*d2c + BETA*ED2V*dk + BETA*EDV*d2k


end subroutine vfuncsp


end module mod_planner

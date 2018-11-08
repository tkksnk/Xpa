module mod_inner


use mod_functions
use mod_spline
implicit none


contains


subroutine inner(mpmat,knotsk,knotsm,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,iterout,vmat0,gmat0,eptime)


    use mod_parameters
    real(8), intent(in) :: mpmat(:,:), knotsk(:), knotsm(:), invTk(:,:), invTm(:,:), Gz(:), Gl(:), Pz(:,:), Gez(:,:), Gy(:), Gd(:), Pez(:,:,:), Py(:,:), Pd(:,:)
    integer, intent(in) :: idmat(:,:), iterout
    real(8), intent(inout) :: vmat0(:,:,:,:), gmat0(:,:,:,:)
    real(8), intent(out) :: eptime
    real(8) znow, lnow, tau, enow, ynow, beta, mnow, know, rr, wr, yterm, klow, khigh, kp, f0, df0, d2f0
    real(8) vcond(nk,nm), vm(nk), cf(4,rk+1), gmat1(nk,nm,nz,nx), vmat1(nk,nm,nz,nx), mp, wm, prob
    real(8), allocatable :: cfmat(:,:,:)
    real(8) crit, diff, diffg, diffv
    integer iter, s1
    integer ik, im, jm, ix, jx, ie, iy, id, je, jy, jd, iz, c1, c2, cr, ii


    call system_clock(c1,cr)

    ! initial guess for the value function
	if (iterout<1) then

        !omp parallel do private(iz,ix,im,ik,ie,iy,znow,lnow,mnow,rr,wr,know,enow,ynow,yterm)
        do iz = 1,nz
            do ix = 1,nx
                do im = 1,nm
                    do ik = 1,nk

                        ie = idmat(ix,1)
                        iy = idmat(ix,2)
                        znow = Gz(iz)
                        lnow = Gl(iz)
                        mnow = knotsm(im)
                        rr = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
                        wr = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

                        know = knotsk(ik)
                        enow = Gez(ie,iz)
                        ynow = Gy(iy)
                        yterm = wr*ynow*enow + rr*know/THETA

                        vmat0(ik,im,iz,ix) = log(yterm)
                        gmat0(ik,im,iz,ix) = know

                    end do
                end do
            end do
        end do

    end if

    if (linflag) then
        allocate(cfmat(4,rk+1,1))
    else
        allocate(cfmat(16,rk+1,rm+1))
    end if

    diff = 1d+4
    iter = 0
    s1   = 0

    do while (diff>critin)

        !$omp parallel do private(ix,ie,iy,id,iz,vcond,prob,jx,je,jy,jd,im,znow,lnow,mnow,mp,rr,wr,jm,wm,vm,cfmat)
        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            do iz = 1,nz

                vcond = 0.0d0
                prob = 0.0d0
                do jx = 1,nx

                    je = idmat(jx,1)
                    jy = idmat(jx,2)
                    jd = idmat(jx,3)
                    vcond = vcond + Py(iy,jy)*Pd(id,jd)*(Pez(ie,je,2*(iz-1)+1)*vmat0(:,:,1,jx) + Pez(ie,je,2*(iz-1)+2)*vmat0(:,:,2,jx))
                    prob = prob + Py(iy,jy)*Pd(id,jd)*(Pez(ie,je,2*(iz-1)+1) + Pez(ie,je,2*(iz-1)+2))

                end do

                if (linflag .eqv. .false.) cfmat = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

                do im = 1,nm

                    znow = Gz(iz)
                    lnow = Gl(iz)
                    mnow = knotsm(im)
                    mp = mpmat(im,iz)
                    rr = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
                    wr = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

                    ! interpolate EVmat over m and fit 1-dim spline
                    if (linflag) then
                        jm = gridlookup2(mp,knotsm)
                        wm = (knotsm(jm+1)-mp)/(knotsm(jm+1)-knotsm(jm))
                        vm = wm*vcond(:,jm) + (1.0d0-wm)*vcond(:,jm+1)
                        cfmat(:,:,1) = spfit(invTk,vm,rk,knotsk)
                    end if

                    !$omp parallel do private(ik,know,enow,ynow,beta,yterm,klow,khigh,kp,f0,df0,d2f0)
                    do ik = 1,nk

                        know = knotsk(ik)
                        enow = Gez(ie,iz)
                        ynow = Gy(iy)
                        beta = Gd(id)
                        yterm = wr*ynow*enow + rr*know/THETA

                        if (s1==0) then

                            klow = knotsk(1)
                            khigh = min(yterm,knotsk(nk))
                            if (nraflag) then
                                call nra(gmat0(ik,im,iz,ix),klow,khigh,kp,mp,cfmat,yterm,knotsk,knotsm,beta)
                            else
                                call gss(klow,khigh,kp,mp,cfmat,yterm,knotsk,knotsm,beta)
                            end if
                            call vfuncsp2(kp,mp,cfmat,yterm,knotsk,knotsm,beta,f0,df0,d2f0)
                            gmat1(ik,im,iz,ix) = kp
                            vmat1(ik,im,iz,ix) = -f0

                        else

                            kp = gmat0(ik,im,iz,ix)
                            call vfuncsp2(kp,mp,cfmat,yterm,knotsk,knotsm,beta,f0,df0,d2f0)
                            gmat1(ik,im,iz,ix) = kp
                            vmat1(ik,im,iz,ix) = -f0

                        end if

                    end do ! loop for k

                end do ! loop for m

            end do ! loop for x

        end do ! loop for z

        diffg = maxval(maxval(maxval(maxval(abs(gmat1-gmat0),1),1),1),1)
        diffv = maxval(maxval(maxval(maxval(abs(vmat1-vmat0),1),1),1),1)
        diff  = diffv
        iter  = iter + 1
        ! diagnosis
        if (mod(iter,diagnum)==0) then

            write(*,"('  iteration ', I4, '  ||Tv-v|| = ', F10.5, '  ||Tg-g|| = ', F10.5)") iter, diffv, diffg

        end if

        if (diffg<1d-4 .and. s1==0) then
           s1 = 1
        elseif (s1>0 .and. s1<20) then
           s1 = s1+1
        elseif (s1>=20) then
           s1 = 0
        end if

        gmat0 = gmat1
        vmat0 = vmat1

    end do


    call system_clock(c2, cr)
    eptime = dble(c2-c1)/cr
    ! write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine inner


subroutine nra(xinit,xmin,xmax,x,mp,cfmat,yterm,knotsk,knotsm,beta)


    use mod_parameters, only: critn
	real(8), intent(in) :: xinit, xmin, xmax, mp, cfmat(:,:,:), yterm, knotsk(:), knotsm(:), beta
	real(8), intent(out) :: x
    real(8) FL, DFL, FH, DFH, f0, df0, d2f0, NewtonStep, x0, xlow, xhigh, x1, diff
    integer iter


    call vfuncsp2(xmin,mp,cfmat,yterm,knotsk,knotsm,beta,FL,DFL,D2F0)
    call vfuncsp2(xmax,mp,cfmat,yterm,knotsk,knotsm,beta,FH,DFH,D2F0)

    if (DFL<0.0d0) then

        x0 = xmin

    elseif (DFH>0.0d0) then

        x0 = xmax

    else ! DFL>=0 amd DFH<=0

        x0    = xinit
        xlow  = xmin
        xhigh = xmax

        diff = 1d+4
        iter = 0

        do while (diff>critn)

            call vfuncsp2(x0,mp,cfmat,yterm,knotsk,knotsm,beta,F0,DF0,D2F0)
            NewtonStep = df0/d2f0
            x1 = x0 - NewtonStep

            if (df0>0.0d0) then
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
            ! diagnosis
            ! write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff
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


subroutine gss(xmin,xmax,x,mp,cfmat,yterm,knotsk,knotsm,beta)


    use mod_parameters, only: critg
	real(8), intent(in) :: xmin, xmax, mp, cfmat(:,:,:), yterm, knotsk(:), knotsm(:), beta
	real(8), intent(out) :: x
	real(8) rg, a, b, c, d, z, fc, df0, d2f0, fd, diff
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    call vfuncsp2(c,mp,cfmat,yterm,knotsk,knotsm,beta,fc,DF0,D2F0)
	d = a + (1.0d0-rg)*(b-a)
    call vfuncsp2(d,mp,cfmat,yterm,knotsk,knotsm,beta,fd,DF0,D2F0)

!	crit = 1d-5
    diff = 1d+4
	iter = 0

	do while (diff>critg)

		if (fc>=fd) then

            z = c + (1.0d0-rg)*(b-c)
           	a = c
           	c = d
           	fc = fd
           	d = z
            call vfuncsp2(d,mp,cfmat,yterm,knotsk,knotsm,beta,fd,DF0,D2F0)

        else

           	z = a + rg*(d-a)
           	b = d
           	d = c
           	fd = fc
           	c = z
            call vfuncsp2(c,mp,cfmat,yterm,knotsk,knotsm,beta,fc,DF0,D2F0)

       	end if

        diff = d-c
        iter = iter+1

        ! diagnosis
        ! write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff

	end do

	x = c


end subroutine gss


subroutine vfuncsp2(kp,mp,cfmat,yterm,knotsk,knotsm,beta,F0,DF0,D2F0)


    use mod_parameters, only: rk, rm, SIGMA, THETA, linflag
	real(8), intent(in) :: kp, mp, cfmat(:,:,:), yterm, knotsk(:), knotsm(:), beta
    real(8), intent(out) :: F0, DF0, D2F0
	real(8) cnow, EV, EDV, ED2V


    cnow = yterm - kp
    if (linflag) then
        call speva(reshape(cfmat,(/4,rk+1/)),kp,rk,knotsk,EV,EDV,ED2V)
    else
        call speva2(cfmat,kp,mp,rk,rm,knotsk,knotsm,EV,EDV,ED2V)
    end if

    if (SIGMA==1.0d0) then
        F0 = log(cnow) + BETA*THETA*ev
    else
        F0 = (cnow**(1.0d0-SIGMA)-1.0d0)/(1.0d0-SIGMA) + BETA*THETA*ev
    end if
    F0   = -F0 ! for gss
    DF0  = -cnow**(-SIGMA) + BETA*THETA*edv
    D2F0 = -SIGMA*cnow**(-SIGMA-1.0d0) + BETA*THETA*ed2v


end subroutine vfuncsp2


end module mod_inner

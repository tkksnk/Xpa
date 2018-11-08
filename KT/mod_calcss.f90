module mod_calcss


use mod_functions
use mod_spline
implicit none


contains


subroutine calcss(knotsk,knotsb,invTk,Ge,Pe,mue,vmat0,gmat0,mu0,evec,mpmat,ymat,nmat)


    use mod_parameters
    real(8), intent(in) :: knotsk(:), knotsb(:), invTk(:,:), Ge(:), Pe(:,:), mue(:)
    real(8), intent(out) :: vmat0(:,:), gmat0(:,:), mu0(:,:), evec(:,:), mpmat(:,:), ymat(:,:), nmat(:,:)
    integer iter, ib, ie, je, ct1, ct2, cr
    real(8) znow, wL, wLnew, BL, wH, wHnew, BH, diff, w0, w1, B0, p0
    real(8) gwmat0(ne), gcmat0(nk,ne), cfmate(4,rk+1,ne)
    real(8) mnow, know, vcond(nk), cf(4,rk+1), enow, kwnow, f0, df0, d2f0, e0now, klow, khigh, kcnow, f1, df1, d2f1, e1now, xinow, alpha, mp, yterm, nnow, ynow, inow, cnow
    real(8) mpvec(ne), yvec(ne), cvec(ne), p1, eptime


    call system_clock(ct1,cr)

    znow = 1.0d0
    ! initial distribution
    ! NOTE: 22 Jan 2018 The initial distribution below is much more stable
	mu0 = 1.0d0/dble(nb*ne)

    if (bisectw) then
    	! bisection method
        ! NOTE: initial lower and upper bounds are given arbitrary
    	wL = 0.90d0
        call driver(wL,znow,knotsk,knotsb,invTk,Ge,Pe,0,wLnew,mpmat,ymat,nmat,vmat0,gmat0,gwmat0,gcmat0,mu0)
        BL = wL-wLnew
        print *, BL

    	wH = 1.10d0
        call driver(wH,znow,knotsk,knotsb,invTk,Ge,Pe,0,wHnew,mpmat,ymat,nmat,vmat0,gmat0,gwmat0,gcmat0,mu0)
        BH = wH-wHnew
        print *, BH

    	diff = 1e+4
    	iter = 0

    	do while(diff>critbp)

    	    w0 = (wL+wH)/2.0d0
            call driver(w0,znow,knotsk,knotsb,invTk,Ge,Pe,iter,w1,mpmat,ymat,nmat,vmat0,gmat0,gwmat0,gcmat0,mu0)
    	    B0 = w0-w1

    	    if (B0*BL>0.0d0) then
    	        wL = w0
    	        BL = B0
    	    else
    	        wH = w0
    	    end if

    	    diff = wH-wL
    	    iter = iter + 1

            write(*,"('  bisection ', I3, '  wH-wL = ', F10.5)") iter, diff
            print *, B0

    	end do

    else
    ! iterative method
        w0 = 1.0d0
        diff = 1e+4
        iter = 0

        do while (diff>critbp)

            call driver(w0,znow,knotsk,knotsb,invTk,Ge,Pe,iter,w1,mpmat,ymat,nmat,vmat0,gmat0,gwmat0,gcmat0,mu0)

            diff = abs(log(w1)-log(w0))
            iter = iter + 1

            write(*,"('  iteration ', I3, '  w1-w0 = ', F10.5)") iter, diff

            ! updating price
            w0 = dampss*w1 + (1.0d0-dampss)*w0

        end do

    end if


    ! calculate the bias correction terms
    p0 = ETA/w0
    mnow = 0.0d0
    do ie = 1,ne

        mnow = mnow + sum(knotsb*mu0(:,ie),1)/sum(mu0(:,ie),1)

        vcond = 0.0d0
        do je = 1,ne

            vcond = vcond + Pe(ie,je)*reshape(vmat0(:,je),(/nk/))

        end do

        cfmate(:,:,ie) = spfit(invTk,vcond,rk,knotsk)

    end do

    do ie = 1,ne

        ! explicit aggregation: aggregate capital indexed by individual productivity
        if (naiveflag) then
            know = mnow
        else
            know = sum(knotsb*mu0(:,ie),1)/sum(mu0(:,ie),1)
        end if

        enow = Ge(ie)
        cf = cfmate(:,:,ie)

        call gss(knotsk(1),knotsk(nk),kwnow,cf,knotsk,p0)
        call vfuncsp(kwnow,cf,knotsk,p0,f0,df0,d2f0)
        e0now = -f0

        klow  = (1.0d0-DELTA-B)/GAMY*know
        khigh = (1.0d0-DELTA+B)/GAMY*know
        call gss(klow,khigh,kcnow,cf,knotsk,p0)
        call vfuncsp(kcnow,cf,knotsk,p0,f1,df1,d2f1)
        e1now = -f1

        xinow = min(XIBAR,max(0.0d0,(e0now-e1now)/ETA))
        alpha = xinow/XIBAR

        mp = alpha*kwnow + (1.0d0-alpha)*kcnow
        yterm = enow*know**THETA
        nnow = (NU*yterm/w0)**(1.0d0/(1.0d0-NU))
        ynow = yterm*nnow**NU
        inow = GAMY*(alpha*kwnow + (1.0d0-alpha)*kcnow) - (1.0d0-DELTA)*know
        cnow = ynow - inow

        mpvec(ie) = mp
        yvec(ie) = ynow
        cvec(ie) = cnow

    end do

    ! calculate the Jensen's inequality
    do ie = 1,ne

        evec(ie,1) = -(mpvec(ie) - sum(mpmat(:,ie)*mu0(:,ie),1)/sum(mu0(:,ie),1))
        evec(ie,2) = -(yvec(ie) - sum(ymat(:,ie)*mu0(:,ie),1)/sum(mu0(:,ie),1))

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine calcss


subroutine driver(w0,znow,knotsk,knotsb,invTk,Ge,Pe,iterout,w1,mpmat,ymat,nmat,vmat0,gmat0,gwmat0,gcmat0,mu0)


    use mod_parameters
    implicit none
    real(8), intent(in) :: w0, znow, knotsk(:), knotsb(:), invTk(:,:), Ge(:), Pe(:,:)
    integer, intent(in) :: iterout
    real(8), intent(inout) :: vmat0(:,:), mu0(:,:)
    real(8), intent(out) :: w1, mpmat(:,:), ymat(:,:), nmat(:,:), gmat0(:,:), gwmat0(:), gcmat0(:,:)
    real(8) pimat(nk,ne), vmat1(nk,ne), gmat1(nk,ne), gwmat1(ne), gcmat1(nk,ne)
    real(8) vcond(nk), cf(4,rk+1), cfmate(4,rk+1,ne)
    real(8) p0, enow, know, yterm, nnow, ynow, klow, khigh, kwnow, f0, df0, d2f0, e0now, kcnow, f1, df1, d2f1, e1now, xinow, alpha
    real(8) crit, diff, diffgw, diffgc, diffg, diffv
    integer iter, s1, s2, ik, ie, je, ct1, ct2, cr
    integer ib, SOLVE(nb,ne), kb1(ne), kb2(nb,ne)
    real(8) G(nb,ne,3), wb1(ne), a1(nb,ne), wb2(nb,ne), inow, imat(nb,ne), ikmat(nb,ne), mu1(nb,ne)

    real(8) time_begin, time_end, eptime


    call system_clock(ct1,cr)
    ! write(*,"('  INNER LOOP at w = ', F10.5)") w0
    p0 = ETA/w0


    ! initial guess for vmat0
    do ie = 1,ne

        do ik = 1,nk

            enow = Ge(ie)
            know = knotsk(ik)
            yterm = znow*enow*know**THETA
            nnow = (NU*yterm/w0)**(1.0d0/(1.0d0-NU))
            ynow = yterm*nnow**NU
            pimat(ik,ie) = (ynow - w0*nnow + (1.0d0-DELTA)*know)*p0

        end do

    end do

    if (iterout<1) vmat0 = pimat
    ! vmat0 = pimat
    gmat0 = 0.0d0
    gwmat0 = 0.0d0
    gcmat0 = 0.0d0

    vmat1 = 0.0d0
    gmat1 = 0.0d0
    gwmat1 = 0.0d0
    gcmat1 = 0.0d0

    diff = 1d+4
    iter = 0
    s1   = 0
    s2   = 0

    do while (diff>critin)

        !$omp parallel do private(vcond,cf,kwnow,f0,df0,d2f0,e0now)
        do ie = 1,ne

            vcond = 0.0d0
            do je = 1,ne

                vcond = vcond + Pe(ie,je)*reshape(vmat0(:,je),(/nk/))

            end do

            cf = spfit(invTk,vcond,rk,knotsk)

            if (s1==0) then

                if (nraflag) then
                    call nra(1.0d0,knotsk(1),knotsk(nk),kwnow,cf,knotsk,p0)
                else
                    call gss(knotsk(1),knotsk(nk),kwnow,cf,knotsk,p0)
                end if

            else

                kwnow = gwmat0(ie)

            end if

            call vfuncsp(kwnow,cf,knotsk,p0,f0,df0,d2f0)
            e0now = -f0 !-vfuncsp(kwnow,cf,knotsk,p0)
            gwmat1(ie) = kwnow

            !$omp parallel do private(know,klow,khigh,kcnow,f1,df1,d2f1,e1now,xinow,alpha)
            do ik = 1,nk

                if (s2==0) then

                    know  = knotsk(ik)
                    klow  = (1.0d0-DELTA-B)/GAMY*know
                    khigh = (1.0d0-DELTA+B)/GAMY*know

                    if (nraflag) then
                        call nra(know,klow,khigh,kcnow,cf,knotsk,p0)
                    else
                        call gss(klow,khigh,kcnow,cf,knotsk,p0)
                    end if

                else

                    kcnow = gcmat0(ik,ie)

                end if

                call vfuncsp(kcnow,cf,knotsk,p0,f1,df1,d2f1)
                e1now = -f1 !-vfuncsp(kcnow,cf,knotsk,p0)
                gcmat1(ik,ie) = kcnow
                xinow = min(XIBAR,max(0.0d0,(e0now-e1now)/ETA))
                alpha = xinow/XIBAR

                vmat1(ik,ie) = pimat(ik,ie) - ETA*xinow**2/2.0d0/XIBAR &
                    + alpha*e0now + (1.0d0-alpha)*e1now
                gmat1(ik,ie) = alpha*kwnow + (1.0d0-alpha)*kcnow

            end do ! loop for k
            !$omp end parallel do

        end do ! loop for e
        !$omp end parallel do

        diffgw = maxval(abs(gwmat1-gwmat0),1)
        diffgc = maxval(maxval(abs(gcmat1-gcmat0),1),1)
		diffg  = max(diffgw,diffgc)
        diffv  = maxval(maxval(abs(vmat1-vmat0),1),1)
        diff   = diffv
        iter = iter + 1

        ! diagnosis
        if (mod(iter,200)==0) then

            write(*,"('  iteration ', I4, '  ||Tv-v|| = ', F10.5, '  ||Tg-g|| = ', F10.5)") iter, diffv, diffg

        end if

        if (diffgw<1d-4 .and. s1==0) then
            s1 = 1
        elseif (s1>0 .and. s1<20) then
            s1 = s1+1
        elseif (s1>=20) then
            s1 = 0
        end if

        if (diffgc<1d-4 .and. s2==0) then
            s2 = 1
        elseif (s2>0 .and. s2<20) then
            s2 = s2+1
        elseif (s2>=20) then
            s2 = 0
        end if

        gwmat0 = gwmat1
        gcmat0 = gcmat1
        vmat0  = vmat1
        gmat0  = gmat1

    end do


    ! fit splines to conditional expectation
    do ie = 1,ne

        vcond = 0.0d0

        do je = 1,ne

            vcond = vcond + Pe(ie,je)*reshape(vmat0(:,je),(/nk/))

        end do

        cfmate(:,:,ie) = spfit(invTk,vcond,rk,knotsk)

    end do


    diff = 1d+4
    iter = 0

    ! NOTE: need to initialize matrices as some points in the rectangle grid are never visited
    mpmat = 0.0d0
    nmat = 0.0d0
    ymat = 0.0d0
    imat = 0.0d0
    ikmat = 0.0d0

    do while(diff>critmu .and. iter<1000)

        ! NOTE: it should be before the while loop?
        SOLVE = 0
        G = 0.0d0

        !$omp parallel do private(enow,cf,kwnow,f0,df0,d2f0,e0now)
        do ie = 1,ne

            enow = Ge(ie)
            cf = reshape(cfmate(:,:,ie), (/4,rk+1/))

            ! NOTE: this part also can be economized
            if (nraflag) then
                call nra(1.0d0,knotsk(1),knotsk(nk),kwnow,cf,knotsk,p0)
            else
                call gss(knotsk(1),knotsk(nk),kwnow,cf,knotsk,p0)
            end if

            call vfuncsp(kwnow,cf,knotsk,p0,f0,df0,d2f0)
            e0now = -f0 !-vfuncsp(kwnow,cf,knotsk,p0)

            ! k1, w1 are index and weight at each grid for kw
            kb1(ie) = gridlookup2(kwnow,knotsb)
            wb1(ie) = (knotsb(kb1(ie)+1)-kwnow)/(knotsb(kb1(ie)+1)-knotsb(kb1(ie)))

            !$omp parallel do private(know,klow,khigh,kcnow,f1,df1,d2f1,e1now,xinow,yterm,nnow,ynow,inow)
            do ib = 1,nb

                know = knotsb(ib)

                ! if the distribution is greater than 0
                if (mu0(ib,ie)>0.0d0) then

                    ! if the grid point is not visited before
                    if (SOLVE(ib,ie)==0) then

                        ! solve for xi
                        klow  = (1.0d0-DELTA-B)/GAMY*know
                        khigh = (1.0d0-DELTA+B)/GAMY*know
                        if (nraflag) then
                            call nra(know,klow,khigh,kcnow,cf,knotsk,p0)
                        else
                            call gss(klow,khigh,kcnow,cf,knotsk,p0)
                        end if

                        call vfuncsp(kcnow,cf,knotsk,p0,f1,df1,d2f1)
                        e1now = -f1 !-vfuncsp(kcnow,cf,knotsk,p0)
                        xinow = min(XIBAR,max(0.0d0,(e0now-e1now)/ETA))
                        ! adjusting probability
                        a1(ib,ie) = xinow/XIBAR
                        ! store the results for later use
                        G(ib,ie,1) = kwnow
                        G(ib,ie,2) = kcnow
                        G(ib,ie,3) = a1(ib,ie)
                        SOLVE(ib,ie) = 1

                    else

                        kwnow = G(ib,ie,1)
                        kcnow = G(ib,ie,2)
                        a1(ib,ie) = G(ib,ie,3)
                        xinow = a1(ib,ie)*XIBAR

                    end if

                    yterm = znow*enow*know**THETA
                    nnow = (NU*yterm/w0)**(1.0d0/(1.0d0-NU))
                    ynow = yterm*nnow**NU
                    inow = GAMY*(a1(ib,ie)*kwnow + (1.0d0-a1(ib,ie))*kcnow) - (1.0d0-DELTA)*know

                    ! distribution
                    mpmat(ib,ie) = a1(ib,ie)*kwnow + (1.0d0-a1(ib,ie))*kcnow
                    nmat(ib,ie) = nnow + xinow**2/XIBAR/2.0d0 ! NOTE: Is the adj. cost term correct?
                    ymat(ib,ie) = ynow
                    imat(ib,ie) = inow
                    ikmat(ib,ie) = inow/know

                    kb2(ib,ie) = gridlookup2(kcnow,knotsb)
                    wb2(ib,ie) = (knotsb(kb2(ib,ie)+1)-kcnow)/(knotsb(kb2(ib,ie)+1)-knotsb(kb2(ib,ie)))

                end if

            end do
            !$omp end parallel do

        end do
        !$omp end parallel do

        ! update the disribution using k1, w1, k2, w2, and a1
        mu1 = 0.0d0
        do ie = 1,ne

            do ib = 1,nb

                if (mu0(ib,ie)>0.0d0) then

                    do je = 1,ne

                        ! if adjust with prob a1(ib,ie)
                        mu1(kb1(ie),je)      = mu1(kb1(ie),je)      + mu0(ib,ie)*wb1(ie)*a1(ib,ie)*Pe(ie,je)
                        mu1(kb1(ie)+1,je)    = mu1(kb1(ie)+1,je)    + mu0(ib,ie)*(1.0d0-wb1(ie))*a1(ib,ie)*Pe(ie,je)
                        ! if not adjust with prob 1-a1(ib,ie)
                        mu1(kb2(ib,ie),je)   = mu1(kb2(ib,ie),je)   + mu0(ib,ie)*wb2(ib,ie)*(1.0d0-a1(ib,ie))*Pe(ie,je)
                        mu1(kb2(ib,ie)+1,je) = mu1(kb2(ib,ie)+1,je) + mu0(ib,ie)*(1.0d0-wb2(ib,ie))*(1.0d0-a1(ib,ie))*Pe(ie,je)

                    end do

                end if

            end do

        end do

        diff = maxval(maxval(abs(mu1-mu0),1),1)
        iter = iter + 1

        ! diagnosis
        if (mod(iter,200)==0) then

            write(*,"('  iteration ', I6, '  ||mu1-mu0|| = ', F15.10)") iter, diff

        end if

        ! update distribution
        mu0 = mu1

    end do

    ! new w
    w1 = ETA*sum(sum(mu0*(ymat-imat),1),1)


    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    ! write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine driver


subroutine gss(xmin,xmax,x,cf,knotsk,p)


    use mod_parameters, only: rk, critg
	real(8), intent(in) :: xmin, xmax, p, knotsk(:), cf(:,:)
	real(8), intent(out) :: x

	real(8) rg, a, b, c, d, z, fc, fd, df, d2f, crit, diff
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    call vfuncsp(c,cf,knotsk,p,fc,df,d2f)
	d = a + (1.0d0-rg)*(b-a)
    call vfuncsp(d,cf,knotsk,p,fd,df,d2f)

    diff = 1d+4
	iter = 0

	do while (diff>critg)

 		if (fc>=fd) then

        	z = c + (1.0d0-rg)*(b-c)
        	a = c
        	c = d
        	fc = fd
        	d = z
            call vfuncsp(d,cf,knotsk,p,fd,df,d2f)

	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            call vfuncsp(c,cf,knotsk,p,fc,df,d2f)

    	end if

	    diff = d-c
    	iter = iter+1

        ! diagnosis
!        write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff

	end do

	x = c


end subroutine gss


subroutine nra(xinit,xmin,xmax,x,cf,knotsk,p)


    use mod_parameters, only: rk, critn
	real(8), intent(in) :: xinit, xmin, xmax, cf(:,:), knotsk(:), p
	real(8), intent(out) :: x
    real(8) FL, DFL, FH, DFH, f0, df0, d2f0, NewtonStep, x0, xlow, xhigh, x1, diff
    integer iter


    call vfuncsp(xmin,cf,knotsk,p,FL,DFL,d2f0)
    call vfuncsp(xmax,cf,knotsk,p,FH,DFH,d2f0)

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

!        Newton-Rhapson
        do while (diff>critn)

            call vfuncsp(x0,cf,knotsk,p,f0,df0,d2f0)
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


subroutine vfuncsp(kp,cf,knotsk,p,f,df,d2f)


    use mod_parameters, only: rk, GAMY, BETA
	real(8), intent(in) :: kp, p, cf(:,:), knotsk(:)
    real(8), intent(out) :: f, df, d2f
    real(8) EV, EDV, ED2V


    call speva(cf,kp,rk,knotsk,EV,EDV,ED2V)
    f = -p*GAMY*kp + BETA*EV
    f = -f
    df  = -p*GAMY + BETA*EDV
    d2f = BETA*ED2V


end subroutine vfuncsp


subroutine calcdiststat(knotsb,mu0,mpmat,ikvec)


    use mod_parameters, only: ne, nb, GAMY, DELTA
    real(8), intent(in) :: knotsb(:), mu0(:,:), mpmat(:,:)
    real(8), intent(out) :: ikvec(:)
    real(8) ikmat(nb,ne)
    integer ie, ib


    do ie = 1,ne

        do ib = 1,nb

            ikmat(ib,ie) = GAMY*mpmat(ib,ie)/knotsb(ib) - (1.0d0-DELTA)

        end do

    end do


    ikvec = 0.0d0 ! NOTE: need to initialize ikvec; otherwise it has values in the previous loop?
    ikvec(1) = sum(sum(mu0*ikmat,1),1)

    do ie = 1,ne

        do ib = 1,nb

            ikvec(2) = ikvec(2) + mu0(ib,ie)*(ikmat(ib,ie)-ikvec(1))**2
            if (abs(ikmat(ib,ie))<0.01d0) ikvec(3) = ikvec(3) + mu0(ib,ie)
            if (ikmat(ib,ie)>0.20d0)      ikvec(4) = ikvec(4) + mu0(ib,ie)
            if (ikmat(ib,ie)<-0.20d0)     ikvec(5) = ikvec(5) + mu0(ib,ie)
            if (ikmat(ib,ie)>=0.01d0)     ikvec(6) = ikvec(6) + mu0(ib,ie)
            if (ikmat(ib,ie)<=-0.01d0)    ikvec(7) = ikvec(7) + mu0(ib,ie)

        end do

    end do


end subroutine calcdiststat


end module mod_calcss

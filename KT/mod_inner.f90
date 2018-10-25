module mod_inner


#ifdef MATLAB_MEX_FILE
USE INTERRUPTEXECUTION
#endif
use mod_functions
use mod_spline
implicit none


contains


subroutine inner(mpmat,pmat,knotsk,knotsm,invTk,invTm,Gz,Pz,Ge,Pe,vmat0,gmat0,eptime,error)


    use mod_parameters
    real(8), intent(in)  :: mpmat(:,:), pmat(:,:), knotsk(:), knotsm(:), invTk(:,:), invTm(:,:), Gz(:), Pz(:,:), Ge(:), Pe(:,:)
    integer, intent(out) :: error
    real(8), intent(out) :: vmat0(:,:,:,:), gmat0(:,:,:,:), eptime
    real(8) pimat(nk,nm,nz,ne), gwmat0(nm,nz,ne), gcmat0(nk,nm,nz,ne), vmat1(nk,nm,nz,ne), gmat1(nk,nm,nz,ne), gwmat1(nm,nz,ne), gcmat1(nk,nm,nz,ne)
    real(8) vcond(nk,nm), cfmat(16,rk+1,rm+1)
    real(8) znow, mnow, mp, p0, w0, enow, know, yterm, nnow, ynow, klow, khigh, kwnow, f0, df0, d2f0, e0now, kcnow, f1, df1, d2f1, e1now, xinow, alpha
    real(8) crit, diff, diffgw, diffgc, diffg, diffv
    integer iter, s1, s2, ik, im, iz, jz, ie, je, ct1, ct2, cr

    real(8) time_begin, time_end

#ifdef MATLAB_MEX_FILE
    integer, external :: mexPrintf
    integer k
    character(len=240) line
    LOGICAL :: INTERRUPTED
#endif


    call system_clock(ct1,cr)

    !$omp parallel do private(znow,mnow,mp,p0,w0,enow,know,yterm,nnow,ynow)
    do ie = 1,ne

        do iz = 1,nz

            znow = Gz(iz)

            do im = 1,nm

                mnow = knotsm(im)
                mp = mpmat(im,iz)
                p0 = pmat(im,iz)
                w0 = ETA/p0

                do ik = 1,nk

                    enow = Ge(ie)
                    know = knotsk(ik)
                    yterm = znow*enow*know**THETA
                    nnow = (NU*yterm/w0)**(1.0d0/(1.0d0-NU))
                    ynow = yterm*nnow**NU
                    pimat(ik,im,iz,ie) = (ynow - w0*nnow + (1.0d0-DELTA)*know)*p0

                end do

            end do

        end do

    end do

    vmat0 = pimat
    gmat0 = 0.0d0
    gwmat0 = 0.0d0
    gcmat0 = 0.0d0

    vmat1 = 0.0d0
    gmat1 = 0.0d0
    gwmat1 = 0.0d0
    gcmat1 = 0.0d0

    diff  = 1d+4
    iter  = 0
    s1    = 0
    s2    = 0
    ! error = 0

    do while (diff>critin .and. iter<1000)

        !$omp parallel do private(vcond,cfmat,mnow,mp,p0,w0,kwnow,f0,df0,d2f0,e0now)
        do ie = 1,ne

            do iz = 1,nz

                vcond = 0.0d0
                do jz = 1,nz
                    do je = 1,ne
                        vcond = vcond + Pz(iz,jz)*Pe(ie,je)*reshape(vmat0(:,:,jz,je),(/nk,nm/))
                    end do
                end do

                cfmat = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

                do im = 1,nm

                    mp = mpmat(im,iz)
                    p0 = pmat(im,iz)
                    w0 = ETA/p0

                    if (s1==0) then

                        if (nraflag) then
                            call nra(knotsm(im),knotsk(1),knotsk(nk),kwnow,cfmat,knotsk,knotsm,mp,p0)
                        else
                            call gss(knotsk(1),knotsk(nk),kwnow,cfmat,knotsk,knotsm,mp,p0)
                        end if

                    else

                        kwnow = gwmat0(im,iz,ie)

                    end if

                    call vfuncsp2(kwnow,cfmat,knotsk,knotsm,mp,p0,f0,df0,d2f0)
                    e0now = -f0 !-vfuncsp2(kwnow,cfmat,knotsk,knotsm,mp,p0)
                    gwmat1(im,iz,ie) = kwnow

                    !$omp parallel do private(know,klow,khigh,kcnow,f1,df1,d2f1,e1now,xinow,alpha)
                    do ik = 1,nk

                        if (s2==0) then

                            know  = knotsk(ik)

                            if (B==0.0d0) then ! no adjustment when B=0

                                kcnow = (1.0d0-DELTA)/GAMY*know

                            else

                                klow  = (1.0d0-DELTA-B)/GAMY*know
                                khigh = (1.0d0-DELTA+B)/GAMY*know
                                if (nraflag) then
                                    call nra(know,klow,khigh,kcnow,cfmat,knotsk,knotsm,mp,p0)
                                else
                                    call gss(klow,khigh,kcnow,cfmat,knotsk,knotsm,mp,p0)
                                end if

                            end if

                        else

                            kcnow = gcmat0(ik,im,iz,ie)

                        end if

                        call vfuncsp2(kcnow,cfmat,knotsk,knotsm,mp,p0,f1,df1,d2f1)
                        e1now = -f1 !-vfuncsp2(kcnow,cfmat,knotsk,knotsm,mp,p0)
                        gcmat1(ik,im,iz,ie) = kcnow
                        xinow = min(XIBAR,max(0.0d0,(e0now-e1now)/ETA))
                        alpha = xinow/XIBAR

                        vmat1(ik,im,iz,ie) = pimat(ik,im,iz,ie) - ETA*xinow**2/2.0d0/XIBAR &
                            + alpha*e0now + (1.0d0-alpha)*e1now
                        gmat1(ik,im,iz,ie) = alpha*kwnow + (1.0d0-alpha)*kcnow

                    end do ! loop for k
                    !$omp end parallel do

                end do ! loop for m

            end do ! loop for z

        end do ! loop for e
        !$omp end parallel do

        diffgw = maxval(maxval(maxval(abs(gwmat1-gwmat0),1),1),1)
        diffgc = maxval(maxval(maxval(maxval(abs(gcmat1-gcmat0),1),1),1),1)
		diffg  = max(diffgw,diffgc)
        diffv  = maxval(maxval(maxval(maxval(abs(vmat1-vmat0),1),1),1),1)
        diff   = diffv
        iter = iter + 1

        ! diagnosis
        if (mod(iter,diagnum)==0) then
#ifdef MATLAB_MEX_FILE
            write(line,"('  iteration ', I4, '  ||Tv-v|| = ', F10.5, '  ||Tkw-kw|| = ', F10.5, '  ||Tkc-kc|| = ', F10.5)") &
            iter, diffv, diffgw, diffgc
    		k = mexPrintf(line//achar(10))
    		call mexEvalString("drawnow")
#else
            write(*,"('  iteration ', I4, '  ||Tv-v|| = ', F10.5, '  ||Tkw-kw|| = ', F10.5, '  ||Tkc-kc|| = ', F10.5)") &
            iter, diffv, diffgw, diffgc
#endif
        end if

#ifdef MATLAB_MEX_FILE
		INTERRUPTED = utIsInterruptPendingInFortran()
		IF (INTERRUPTED) RETURN
#endif

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


    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine inner


subroutine gss(xmin,xmax,x,cfmat,knotsk,knotsm,mp,p0)


    use mod_parameters, only: rk, rm, critg
	real(8), intent(in) :: xmin, xmax, knotsk(:), knotsm(:), cfmat(:,:,:), mp, p0
	real(8), intent(out) :: x
	real(8) rg, a, b, c, d, z, fc, fd, df0, d2f0, crit, diff
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    call vfuncsp2(c,cfmat,knotsk,knotsm,mp,p0,fc,df0,d2f0)
	d = a + (1.0d0-rg)*(b-a)
    call vfuncsp2(d,cfmat,knotsk,knotsm,mp,p0,fd,df0,d2f0)

    diff = 1d+4
	iter = 0

	do while (diff>critg)

 		if (fc>=fd) then

        	z = c + (1-rg)*(b-c)
        	a = c
        	c = d
        	fc = fd
        	d = z
            call vfuncsp2(d,cfmat,knotsk,knotsm,mp,p0,fd,df0,d2f0)
            ! fd = vfuncsp2(d,cmat,knotsk,knotsm,mp,p)

	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            call vfuncsp2(c,cfmat,knotsk,knotsm,mp,p0,fc,df0,d2f0)
            ! fc = vfuncsp2(c,cmat,knotsk,knotsm,mp,p)

    	end if

	    diff = d-c
    	iter = iter+1

        ! diagnosis
!        write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff

	end do

	x = c


end subroutine gss


subroutine nra(xinit,xmin,xmax,x,cfmat,knotsk,knotsm,mp,p0)


    use mod_parameters, only: rk, rm, critn
	real(8), intent (in) :: xinit, xmin, xmax, cfmat(:,:,:), knotsk(:), knotsm(:), mp, p0
	real(8), intent (out) :: x
    real(8) DFL, DFH, f0, df0, d2f0, NewtonStep, x0, xlow, xhigh, x1, diff
    integer iter


    call vfuncsp2(xmin,cfmat,knotsk,knotsm,mp,p0,f0,DFL,d2f0)
    call vfuncsp2(xmax,cfmat,knotsk,knotsm,mp,p0,f0,DFH,d2f0)

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

        ! Newton-Rhapson
        do while (diff>critn)

            call vfuncsp2(x0,cfmat,knotsk,knotsm,mp,p0,f0,df0,d2f0)
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


! function vfuncsp2(kp,cmat,knotsk,knotsm,mp,p) result(f)
!
!
!     use mod_parameters, only: rk, rm, GAMY, BETA
! 	real(8), intent (in) :: kp, mp, p
!     real(8), intent (in) :: cmat(16,rk+1,rm+1)
!     real(8), intent (in) :: knotsk(rk+2), knotsm(rm+2)
!     real(8) f
! 	real(8) EV, dfx, dfy, dfz
!
!
!     call speva2(cmat,kp,mp,rk,rm,knotsk,knotsm,EV,dfx,dfy)
!     f = -p*GAMY*kp + BETA*EV ! EV must be a scalar
!     f = -f
!
!
! end function vfuncsp2


subroutine vfuncsp2(kp,cfmat,knotsk,knotsm,mp,p0,f0,df0,d2f0)


    use mod_parameters, only: rk, rm, GAMY, BETA
	real(8), intent (in) :: kp, cfmat(:,:,:), knotsk(:), knotsm(:), mp, p0
    real(8), intent (out) :: f0, df0, d2f0
	real(8) EV, EDV, ED2V


    call speva2(cfmat,kp,mp,rk,rm,knotsk,knotsm,EV,EDV,ED2V)
    f0   = -p0*GAMY*kp + BETA*EV
    f0   = -f0
    df0  = -p0*GAMY + BETA*EDV
    d2f0 = BETA*ED2V


end subroutine vfuncsp2


end module mod_inner

module mod_inner


use mod_utils
use mod_spline
use mod_policies
use mod_parameters


implicit none


contains


subroutine inner(agrid,mgrid,invTa,invTm,bgrid,xgrid,zgrid,pbij,pxij,pzij,ALOW,Yagg,idmat, &
	kappakp,kappaw,kappar,Vmat,Wmat,Nmat,apmat,awpmat,anpmat,lmat)


	real(8), intent(in) :: agrid(:), mgrid(:), invTa(:,:), invTm(:,:), bgrid(:), xgrid(:), zgrid(:), &
		pbij(:,:), pxij(:,:), pzij(:,:), ALOW, Yagg, kappakp(:,:), kappaw(:,:), kappar(:,:)
	integer, intent(in) :: idmat(:,:)
	real(8), intent(inout) :: Vmat(:,:,:,:)
	real(8), intent(out) :: Wmat(:,:,:,:), Nmat(:,:,:,:)
	real(8), intent(out) :: apmat(:,:,:,:), awpmat(:,:,:,:), anpmat(:,:,:,:), lmat(:,:,:,:)
    real(8), allocatable :: TVmat(:,:,:,:), TWmat(:,:,:,:), TNmat(:,:,:,:)
	real(8), allocatable :: apmatnew(:,:,:,:), awpmatnew(:,:,:,:), anpmatnew(:,:,:,:), lmatnew(:,:,:,:)
    real(8), allocatable :: cfmat(:,:,:), EVmat0(:,:,:,:), EVmat(:,:), EV(:), dist_vecd(:), EVmatz(:,:,:), EVmate(:,:)
	integer ib, jb, ix, jx, ie, je, ia, iz, jz, im, jm, iteration, i, kv, ka, km, ke, kz, dampflag, howard, howardcnt, jghz, nghz
    real(8) znow, betanow, xnow, mnow, mp, w, r, Lagg, anow, awhigh, anhigh, ap, awp, anp, &
		wnmax, wnmin, distancemax, distancemean, distance, atilde, TW, TN, Trw, Trn, distancevec(maxiter), dist_nonzero, wm, damp, &
		distpolmax(3), distpolmean(3), zp, wz
    integer ct1, ct2, cr
    real(8) eptime


	! NOTE: can be declared w/o allocatiable arrays
	allocate(TVmat(na,nm,ne,nz),TWmat(na,nm,ne,nz),TNmat(na,nm,ne,nz),EVmat0(na,nm,ne,nz),EVmat(na,nm),EV(na),dist_vecd(na*nm*ne*nz))
	allocate(EVmatz(na,nm,ne),EVmate(na,nm))
	allocate(apmatnew(na,nm,ne,nz),awpmatnew(na,nm,ne,nz),anpmatnew(na,nm,ne,nz),lmatnew(na,nm,ne,nz))

	if (linflag==0) then
		allocate(cfmat(16,ra+1,rm+1))
	elseif (linflag==2) then ! linear interpolation on m
		allocate(cfmat(4,ra+1,1))
	endif


	print *, "  INNER LOOP"
	call system_clock(ct1, cr)

	! initial guess for the value function
	!if ((vmatini==0) .or. (iterout==1)) then
    	Vmat = 0.0d0
	!end if

	damp = 0.0d0
	dampflag = 0
	howard = 0
	howardcnt = 0


	! INNER LOOP: value function iteration
    distance = 1d+4
	apmat = 0.0d0
	awpmat = 0.0d0
	anpmat = 0.0d0
	lmat = 0.0d0
    ! iteration = 1


    ! do while (distance > tol)
	do iteration = 1,maxiter

        !omp parallel do default(shared) private(ie,ib,ix,iz,EVmat,jz,je,jb,jx)
		! compute conditional expected value functions on grid points
		do ie = 1,ne

            ib = idmat(ie,1)
            ix = idmat(ie,2)

			do iz = 1,nz

				! if (conzflag==2) then ! continuous z
				!
				! 	znow = zgrid(iz)
				! 	EVmat = 0.0d0
				!
				! 	do jghz = 1,nghz
				!
				! 		zp = exp(RHOZ*log(znow) + xghz(jghz))
				! 		call findjw(zp,zgrid,jz,wz)
				! 		EVmatz = wz*reshape(Vmat(:,:,:,jz),(/na,nm,ne/)) + (1.0d0-wz)*reshape(Vmat(:,:,:,jz+1),(/na,nm,ne/))
				!
				! 		EVmate = 0.0d0
				! 		do je = 1,ne
				!
				! 			jb = idmat(je,1)
				! 			jx = idmat(je,2)
				! 			EVmate = EVmate + pbij(ib,jb)*pxij(ix,jx)*reshape(EVmatz(:,:,je),(/na,nm/))
				!
				! 		end do
				!
				! 		EVmat = EVmat + wghz(jghz)*EVmate
				!
				! 	end do
				!
				! else

					EVmat = 0.0d0

					do jz = 1,nz

						do je = 1,ne

							jb = idmat(je,1)
							jx = idmat(je,2)
							EVmat = EVmat + pbij(ib,jb)*pxij(ix,jx)*pzij(iz,jz)*reshape(Vmat(:,:,je,jz),(/na,nm/))

						end do

					end do

				! end if

				EVmat0(:,:,ie,iz) = EVmat

            end do

        end do
        !omp end parallel do


        !$omp parallel do default(shared) private(ie,iz,znow,ib,ix,betanow,xnow,EVmat,cfmat,jm,wm,EV,im,mnow,mp,w,r,atilde,ia,anow,awp,anp,TW,TN)
		do ie = 1,ne

            ib = idmat(ie,1)
			ix = idmat(ie,2)
			betanow = bgrid(ib)
			xnow = xgrid(ix)

            do iz = 1,nz

                znow = zgrid(iz)

				EVmat = EVmat0(:,:,ie,iz)
				! fit 2-dim spline
				if (linflag==0) cfmat = spfit2(invTa,invTm,EVmat,ra,rm,agrid,mgrid)

                do im = 1,nm

                    mnow = mgrid(im)

                    ! if (conzflag==1) then

                        mp = exp(kappakp(1,iz) + kappakp(2,iz)*log(mnow))
						! print *, mp
					    ! interpolate EVmat over m and fit 1-dim spline
					    if (linflag==2) then
						    call findiw(mp,mgrid,jm,wm)
						    EV = wm*EVmat(:,jm) + (1.0d0-wm)*EVmat(:,jm+1)
						    cfmat(:,:,1) = spfit(invTa,EV,ra,agrid)
					    end if

						! NOTE: 031818 indices in kappas are changed (iz is now in columns)
					    if (fcsteqn==1) then
						    w = exp(kappaw(1,iz) + kappaw(2,iz)*log(mnow))
                    	    r = (znow**(1.0d0/ALPHA))*ALPHA*((w/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
					    elseif (fcsteqn==2) then
						    r = exp(kappar(1,iz) + kappar(2,iz)*log(mnow)) - DELTA
						    w = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((r+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
					    else
						    w = exp(kappaw(1,iz) + kappaw(2,iz)*log(mnow))
						    r = exp(kappar(1,iz) + kappar(2,iz)*log(mnow)) - DELTA
                        end if

                    ! elseif (conzflag==3) then
					!
                    !     mp = exp(kappakp(1,1) + kappakp(2,1)*log(mnow) + kappakp(3,1)*log(znow) + kappakp(4,1)*log(znow)*log(mnow))
                    !     ! interpolate EVmat over m and fit 1-dim spline
					!     if (linflag==2) then
					! 	    call findiw(mp,mgrid,jm,wm)
					! 	    EV = wm*EVmat(:,jm) + (1.0d0-wm)*EVmat(:,jm+1)
					! 	    cfmat(:,:,1) = spfit(invTa,EV,ra,agrid)
					!     end if
					!
					!     if (fcsteqn==1) then
					! 	    w = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow) + kappaw(4,1)*log(znow)*log(mnow))
                    ! 	    r = (znow**(1.0d0/ALPHA))*ALPHA*((w/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
					!     elseif (fcsteqn==2) then
					! 	    r = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow) + kappar(4,1)*log(znow)*log(mnow)) - DELTA
					! 	    w = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((r+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
					!     else
					! 	    w = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow) + kappaw(4,1)*log(znow)*log(mnow))
					! 	    r = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow) + kappar(4,1)*log(znow)*log(mnow)) - DELTA
                    !     end if
					!
                    ! else ! conzflag = 0 or 2
					!
                    !     mp = exp(kappakp(1,1) + kappakp(2,1)*log(mnow) + kappakp(3,1)*log(znow))
                    !     ! interpolate EVmat over m and fit 1-dim spline
					!     if (linflag==2) then
					! 	    call findiw(mp,mgrid,jm,wm)
					! 	    EV = wm*EVmat(:,jm) + (1.0d0-wm)*EVmat(:,jm+1)
					! 	    cfmat(:,:,1) = spfit(invTa,EV,ra,agrid)
					!     end if
					!
					!     if (fcsteqn==1) then
					! 	    w = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow))
                    ! 	    r = (znow**(1.0d0/ALPHA))*ALPHA*((w/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
					!     elseif (fcsteqn==2) then
					! 	    r = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow)) - DELTA
					! 	    w = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((r+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
					!     else
					! 	    w = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow))
					! 	    r = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow)) - DELTA
                    !     end if
					!
                    ! end if

                    if (Trscale == 0.0d0) then
                        ! for a below atilde, not-working is not feasible (consumption negative)
                        atilde = (alow)/(1.0d0+r)
                    else
						! NOTE: Steady state output is used (it's okay here) WHY???
                        call transfers(0.0d0,0.0d0,r,w,0.0d0,Yagg,Trscale,Transprog,T0,Trn)
                        atilde = (alow - Trn)/(1.0d0+r)
                    end if

                    do ia = 1,na

                        anow = agrid(ia)
						awp = awpmat(ia,im,ie,iz)
						anp = anpmat(ia,im,ie,iz)

						call calcwnval(xnow,anow,r,w,mp,ALOW,Yagg,atilde,betanow,agrid,mgrid,cfmat,EVmat,awp,anp,TW,TN,howard)

						TWmat(ia,im,ie,iz) = TW
						TNmat(ia,im,ie,iz) = TN
						awpmatnew(ia,im,ie,iz) = awp
						anpmatnew(ia,im,ie,iz) = anp

                        if (TWmat(ia,im,ie,iz)>TNmat(ia,im,ie,iz)) then
                            TVmat(ia,im,ie,iz) = TW
                            apmatnew(ia,im,ie,iz) = awp
                            lmatnew(ia,im,ie,iz) = 1
                        else
                            TVmat(ia,im,ie,iz) = TN
                            apmatnew(ia,im,ie,iz) = anp
                            lmatnew(ia,im,ie,iz) = 0
                        end if

                    end do

                end do

            end do

        end do
		!omp end parallel do

        ! wnmax = maxval(maxval(maxval(maxval(TWmat-TNmat,1),1),1),1)
        ! wnmin = minval(minval(minval(minval(TWmat-TNmat,1),1),1),1)

    	distancemax = maxval(maxval(maxval(maxval(abs(TVmat-Vmat),1),1),1),1)
		distancemean = sum(sum(sum(sum(abs(TVmat-Vmat),1),1),1),1)
		!distpolmax(1) = maxval(maxval(maxval(maxval(abs(apmatnew-apmat),1),1),1),1)
		!distpolmean(1) = sum(sum(sum(sum(abs(apmatnew-apmat),1),1),1),1)
		distpolmax(1) = maxval(maxval(maxval(maxval(abs(awpmatnew-awpmat),1),1),1),1)
		distpolmean(1) = sum(sum(sum(sum(abs(awpmatnew-awpmat),1),1),1),1)
		distpolmax(2) = maxval(maxval(maxval(maxval(abs(anpmatnew-anpmat),1),1),1),1)
		distpolmean(2) = sum(sum(sum(sum(abs(anpmatnew-anpmat),1),1),1),1)
		distpolmax(3) = maxval(maxval(maxval(maxval(abs(lmatnew-lmat),1),1),1),1)
		distpolmean(3) = sum(sum(sum(sum(abs(lmatnew-lmat),1),1),1),1)
		dist_vecd = reshape(abs(TVmat-Vmat),(/na*nm*ne*nz/))
	    kv = maxloc(dist_vecd,1)
	    ka = mod(kv-1,na) + 1
	    km = mod((kv-ka)/na,nm) + 1
	    ke = mod(((kv-ka)/na-km+1)/nm,ne) + 1
	    kz = (((kv-ka)/na-km+1)/nm-ke+1)/ne + 1
		distance = TVmat(ka,km,ke,kz)-Vmat(ka,km,ke,kz)
		dist_nonzero = dble(count(dist_vecd>tol)) !/dble(na*nm*ne*nz)

		distancevec(iteration) = distance ! sum (not max) decreases monotonically

        Vmat = damp*Vmat + (1.0d0-damp)*TVmat
		Wmat = damp*Wmat + (1.0d0-damp)*TWmat
		Nmat = damp*Nmat + (1.0d0-damp)*TNmat
		apmat = apmatnew
		awpmat = awpmatnew
		anpmat = anpmatnew
		lmat = lmatnew

        ! diagnosis
        if (mod(iteration,1)==0) then
            ! write(*,"('  iteration ', I4, '   ||TV-V|| = ', F8.5, ' max(W-N) = ', F8.5, ' min(W-N) = ', F8.5)") &
            ! iteration, distance, wnmax, wnmin
			write(*,"('  iteration ', I4, '   ||TV-V|| = ', F10.8, '   sum|TV-V| = ', F10.8)") iteration, distance, distancemean
			! write(*,"('          at ( ', I4, ', ', I2, ', ', I2, ', ', I2, ') ', I6)") ka, km, ke, kz, int(dist_nonzero)
			! !write(*,"('          ||Tg-g||  = ( ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ')')") distpolmax
			! !write(*,"('          sum|Tg-g| = ( ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ')')") distpolmean
			! write(*,"('          ||Tawp-awp||   ||Tanp-anp||   ||Tlmat-lmat|| = ( ', F8.5, ', ', F8.5, ', ', F8.5, ')')") distpolmax
			! write(*,"('          sum|Tawp-awp|  sum|Tanp-anp|  sum|Tlmat-lmat|= ( ', F8.5, ', ', F8.5, ', ', F8.5, ')')") distpolmean
			! write(*,"('          Howard''s improvement algorithm: ', I2)") howardcnt
			! print *, howard
        end if

		!if ((iteration>1) .and. (distancevec(iteration-1)<distancevec(iteration))) print *, " WARNING: Value function is NOT monotonically decreasing"

        if ((distancemax<=tol) .or. (iteration==maxiter)) then

			print *, " "
			if (iteration<maxiter) then
				print *, "  Value function iteration converged:"
			else
				print *, "  Maximum number of iterations reached:"
			end if
			do i=max(iteration-4,1),iteration
				write(*,"('  iteration ', I4, '   ||TV-V|| = ', F10.8)") &
            	i, distancevec(i)
			end do

			exit

		end if

		if (howflag==1) then

			!if ((maxval(distpolmax,1)<=5d-2) .and. (howard==0)) then
			if (howard==0) then
				! print *, "  Howard''s improvement algorithm starts"
				! pause
				howard = 1
			end if

			if (howard==1) then

				howardcnt = howardcnt + 1
				if (howardcnt>iteration/10) then
					howard = 0
					howardcnt = 0
				end if

			end if

		end if

		! if ((distancemax<=1d-2) .and. (dampflag==0)) then
		!
		! 	print *, "  dampling starts..."
		! 	pause
		! 	! zgrid = zgrid0
		! 	! pzij = pzij0
		! 	damp = 0.5d0
		! 	dampflag = 1
		!
		! end if

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Innerloop elasped time = ', F10.5)") eptime


end subroutine inner


subroutine calcwnval(xnow,anow,r,w,mp,ALOW,Yagg,atilde,betanow,agrid,mgrid,cfmat,EVmat,awp,anp,TW,TN,howard)


	real(8), intent(in) :: xnow, anow, r, w, mp, ALOW, Yagg, atilde, betanow
	real(8), intent(in) :: agrid(:), mgrid(:), cfmat(:,:,:), EVmat(:,:)
	integer, intent(in) :: howard
	real(8), intent(inout) :: awp, anp
	real(8), intent(out) :: TW, TN
	real(8) Trw, awhigh, Trn, anhigh


	! value of work
	call transfers(xnow,anow,r,w,HBAR,Yagg,Trscale,Transprog,T0,Trw)
	awhigh = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Trw
	if ((1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow - alow <= 0.0d0) then
		write(*,"('  max consumption (when work) is negative', F10.2)") awhigh
	end if

	if (howard==0) call gss(ALOW,awhigh,awp,tol**2,1,linflag,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Trw,betanow,B0,Hbar,taul)
	TW = -objneg_w(awp,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Trw,betanow,B0,Hbar,taul,linflag)

	! value of not working
	if (anow>atilde) then
		call transfers(xnow,anow,r,w,0.0d0,Yagg,Trscale,Transprog,T0,Trn)
		anhigh = (1.0d0+r)*anow + Trn
		if (howard==0) call gss(ALOW,anhigh,anp,tol**2,0,linflag,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Trn,betanow,B0,Hbar,taul)
		TN = -objneg_n(anp,anow,agrid,mgrid,cfmat,EVmat,mp,r,Trn,betanow,linflag)
	else
		anp = ALOW
		TN = -1d+20
	end if


end subroutine calcwnval


subroutine gss(xmin,xmax,x,critg,workflag,linflag,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,BETA,B0,HBAR,taul)


	real(8), intent(in) :: xmin, xmax, critg, anow, xnow, mp, w, r, Tr, beta, B0, Hbar, taul
    real(8), intent(in) :: agrid(:), mgrid(:), cfmat(:,:,:), EVmat(:,:)
    integer, intent(in) :: workflag, linflag
	real(8), intent(out) :: x

	real(8) rg, a, b, c, d, z, fc, fd, crit, diff, fa, fb
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    if (workflag==1) then
        fc = objneg_w(c,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,BETA,B0,HBAR,taul,linflag)
        fa = objneg_w(a,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,BETA,B0,HBAR,taul,linflag)
        fb = objneg_w(b,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,BETA,B0,HBAR,taul,linflag)
    else
        fc = objneg_n(c,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,BETA,linflag)
        fa = objneg_n(a,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,BETA,linflag)
        fb = objneg_n(a,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,BETA,linflag)
    end if

	d = a + (1.0d0-rg)*(b-a)
    if (workflag==1) then
        fd = objneg_w(d,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,beta,B0,Hbar,taul,linflag)
    else
        fd = objneg_n(d,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,beta,linflag)
    end if

    diff = 1d+4
	iter = 0

	do while (diff>critg)

 		if (fc>=fd) then

        	z = c + (1.0d0-rg)*(b-c)
        	a = c
        	c = d
        	fc = fd
        	d = z
            if (workflag==1) then
                fd = objneg_w(d,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,beta,B0,Hbar,taul,linflag)
            else
                fd = objneg_n(d,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,beta,linflag)
            end if

	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            if (workflag==1) then
                fc = objneg_w(c,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,beta,B0,Hbar,taul,linflag)
            else
                fc = objneg_n(c,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,beta,linflag)
            end if

		end if

		diff = d-c
		iter = iter+1

        ! diagnosis
!        write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff

	end do


    If (fa==max(fa,fb,fc)) then
        x = a
    else if (fb == max(fa,fb,fc)) then
        x = c
    else
	    x = (c+d)/2.0d0
    end if


end subroutine gss


function objneg_w(ap,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Tr,BETA,B0,HBAR,taul,linflag) result(f)


    real(8), intent(in) :: ap, anow, xnow, mp, w, r, Tr, beta, B0, HBAR, taul
    real(8), intent(in) :: agrid(:), mgrid(:), cfmat(:,:,:), EVmat(:,:)
	integer, intent(in) :: linflag
    real(8) cnow, ev, edva, edvm, f
	integer ra, rm

	ra = size(agrid,1)-2
	rm = size(mgrid,1)-2

    cnow = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Tr - ap
	cnow = max(cnow,1d-4)

	if (linflag==1) then
		call lineva2(agrid,mgrid,EVmat,ap,mp,ev)
	elseif (linflag==2) then
		call speva(cfmat,ap,ra,agrid,ev,edva,edvm)
	else
		call speva2(cfmat,ap,mp,ra,rm,agrid,mgrid,ev,edva,edvm)
	end if
	! if (isnan(ev)) then
	! 	print *, ap, agrid(1), agrid(48), mp, mgrid(1), mgrid(5), ev
	! 	pause
	! end if
    f = log(cnow) - B0 + beta*ev
    f = -f


end function objneg_w


function objneg_n(ap,anow,agrid,mgrid,cfmat,EVmat,mp,r,Tr,beta,linflag) result(f)


    real(8), intent(in) :: ap, anow, mp, r, Tr, beta
    real(8), intent(in) :: agrid(:), mgrid(:), cfmat(:,:,:), EVmat(:,:)
	integer, intent(in) :: linflag
    real(8) cnow, ev, edva, edvm, f
	integer ra, rm

	ra = size(agrid,1)-2
	rm = size(mgrid,1)-2

    cnow = (1.0d0+r)*anow + Tr - ap
	cnow = max(cnow,1d-4)
	if (linflag==1) then
		call lineva2(agrid,mgrid,EVmat,ap,mp,ev)
	elseif (linflag==2) then
		call speva(cfmat,ap,ra,agrid,ev,edva,edvm)
	else
		call speva2(cfmat,ap,mp,ra,rm,agrid,mgrid,ev,edva,edvm)
	end if
	! if (isnan(ev)) then
	! 	print *, ap, agrid(1), agrid(48), mp, mgrid(1), mgrid(5), ev
	! 	pause
	! end if
    f = log(cnow) + beta*ev
    f = -f


end function objneg_n


end module mod_inner

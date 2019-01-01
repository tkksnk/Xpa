module mod_outer


use mod_spline
use mod_inner
use mod_policies
implicit none


contains


! subroutine outer(agrid,kgrid,mgrid,invTa,invTm,bgrid,xgrid,zgrid,pbij,pxij,pzij,simTT,izvec,zvec,idmat,Yss, &
! 	kappakp,kappaw,kappar,xghz,wghz,RHOZ,conzflag,Vmat,ALPHA,DELTA,HBAR,ALOW,B0,BETA,taul,wgtT,a0,a1,a2,Trscale,Transprog1,Transprog2,T0, &
! 	tol,tolmkt,outermkt,foreind,bsctr,bspct,maxcountmkt,calcma,bsfixr,diagnum,fcsteqn,linflag,ymat,mumat)
subroutine outer(agrid,kgrid,mgrid,invTa,invTm,bgrid,xgrid,zgrid,pbij,pxij,pzij,simTT,izvec,ALOW,Yss,idmat, &
	kappakp,kappaw,kappar,Vmat,aggsim,mumat)


	real(8), intent(in) :: agrid(:), kgrid(:), mgrid(:), invTa(:,:), invTm(:,:), bgrid(:), xgrid(:), zgrid(:), &
		pbij(:,:), pxij(:,:), pzij(:,:), ALOW, Yss, kappakp(:,:), kappaw(:,:), kappar(:,:), Vmat(:,:,:,:)
	! real(8), intent(in) :: ALPHA, DELTA, HBAR, ALOW, B0, BETA, taul, RHOZ, wgtT, a0, a1, a2, Trscale, Transprog1,Transprog2, T0, &
	! Yss, tol, tolmkt
	! NOTE: 180305 bspct is updated inside the subroutine; 28 May 2018 it is not updated due to stability
	!real(8), intent(inout) :: bspct
    ! real(8), intent(in) :: bspct
	! NOTE: 180302 tolout, diffout are not needed anymore
	! NOTE: tolout, diffout are needed because market-clearing is imposed if diffout<10*tolout
	integer, intent(in) :: simTT, izvec(:), idmat(:,:) !, outermkt, foreind, bsctr, maxcountmkt, calcma, bsfixr, diagnum, fcsteqn, linflag, conzflag
	real(8), intent(out) :: aggsim(:,:)
	real(8), intent(inout) :: mumat(:,:)
	real(8), allocatable :: Tmumat(:,:), hdist(:,:), cdist(:,:), ldist(:,:), taxdist(:,:), transdist(:,:), minasseta(:), minassetd(:)
    real(8), allocatable :: cfmat0(:,:,:,:,:), EVmat(:,:), EVmat0(:,:,:,:), EVmatz(:,:,:), EVmate(:,:)
    integer ib, jb, ix, jx, ie, je, iz, jz, ik, time, countmkt !, jghz, nghz
	real(8) znow, mnow, mp, wf, rf, w0, r0, Lagg0, Ls, Cagg, Hs, Taxagg, Transagg, Ks, Ys, Xagg, Gagg, wfvec(simTT), w0vec(simTT), zp, wz
	real(8) damp, diffma
	integer ct1, ct2, cr, cto1, cto2, cro
    real(8) eptime


	! NOTE: can be declared w/o allocatiable arrays
	allocate(Tmumat(nk,ne),hdist(nk,ne),cdist(nk,ne),ldist(nk,ne),taxdist(nk,ne),transdist(nk,ne),minasseta(nx),minassetd(nx))
	allocate(EVmatz(na,nm,ne),EVmate(na,nm),EVmat(na,nm))


	print *, "  OUTER LOOP"
	call system_clock(cto1, cro)

	! compute conditional expected value functions on grid points
	if (conzflag==2) then

		! NOTE: conditional expectation is calculated in every time (z is continuous)
		allocate(EVmat0(na,nm,ne,1),cfmat0(16,ra+1,rm+1,ne,1))

	else

		allocate(EVmat0(na,nm,ne,nz),cfmat0(16,ra+1,rm+1,ne,nz))

	    do ie = 1,ne

			ib = idmat(ie,1)
			ix = idmat(ie,2)

	        do iz = 1,nz

				EVmat = 0.0d0
	            do jz = 1,nz

	                do je = 1,ne

						jb = idmat(je,1)
						jx = idmat(je,2)
	                    EVmat = EVmat + pbij(ib,jb)*pxij(ix,jx)*pzij(iz,jz)*reshape(Vmat(:,:,je,jz),(/na,nm/))

	                end do

				end do

				EVmat0(:,:,ie,iz) = EVmat
				if (linflag==0) cfmat0(:,:,:,ie,iz) = spfit2(invTa,invTm,EVmat,ra,rm,agrid,mgrid)

	        end do

	    end do

	end if


	! initial range for bisection (the previous error is used after time 2)
	! bspct = bspct0

    do time = 1,simTT

		call system_clock(ct1, cr)
        mumat = mumat/sum(sum(mumat,1),1)

		! if (conzflag==2) then
		!
		! 	! NOTE: conditional expectation is calculated in every time (z is continuous)
		! 	znow = zvec(time)
        !     iz = 1
		!
		! 	! compute conditional expected value functions on grid points
		!     do ie = 1,ne
		!
		! 		ib = idmat(ie,1)
		! 		ix = idmat(ie,2)
		!
		! 		! nz is the number of gaussian quadrature
		! 		EVmat = 0.0d0
		! 		do jghz = 1,nghz
		!
		! 			zp = exp(RHOZ*log(znow) + xghz(jghz))
		! 			call findjw(zp,zgrid,jz,wz)
		! 			EVmatz = wz*reshape(Vmat(:,:,:,jz),(/na,nm,ne/)) + (1.0d0-wz)*reshape(Vmat(:,:,:,jz+1),(/na,nm,ne/))
		!
		! 			EVmate = 0.0d0
		! 			do je = 1,ne
		!
		! 				jb = idmat(je,1)
		! 				jx = idmat(je,2)
		! 				EVmate = EVmate + pbij(ib,jb)*pxij(ix,jx)*reshape(EVmatz(:,:,je),(/na,nm/))
		!
		! 			end do
		!
		! 			EVmat = EVmat + wghz(jghz)*EVmate
		!
		! 		end do
		!
		! 		EVmat0(:,:,ie,1) = EVmat
		! 		if (linflag==0) cfmat0(:,:,:,ie,1) = spfit2(invTa,invTm,EVmat,ra,rm,agrid,mgrid)
		!
		!     end do
		!
		! else

			iz = izvec(time)
	        znow = zgrid(iz)

		! end if

        mnow = 0.0d0
        do ie = 1,ne
            mnow = mnow + sum(mumat(:,ie)*kgrid)
        end do

        ! if (conzflag==1) then

			mp = exp(kappakp(1,iz) + kappakp(2,iz)*log(mnow))
			! NOTE: 031818 indices in kappas are changed (iz is now in columns)
			if (fcsteqn==1) then
				wf = exp(kappaw(1,iz) + kappaw(2,iz)*log(mnow))
				rf = (znow**(1.0d0/ALPHA))*ALPHA*((wf/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
			elseif (fcsteqn==2) then
				rf = exp(kappar(1,iz) + kappar(2,iz)*log(mnow)) - DELTA
				wf = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rf+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
			else
				wf = exp(kappaw(1,iz) + kappaw(2,iz)*log(mnow))
				rf = exp(kappar(1,iz) + kappar(2,iz)*log(mnow)) - DELTA
			end if

        ! elseif (conzflag==3) then
		!
		! 	mp = exp(kappakp(1,1) + kappakp(2,1)*log(mnow) + kappakp(3,1)*log(znow) + kappakp(4,1)*log(znow)*log(mnow))
		! 	if (fcsteqn==1) then
		! 		wf = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow) + kappaw(4,1)*log(znow)*log(mnow))
		! 		rf = (znow**(1.0d0/ALPHA))*ALPHA*((wf/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
		! 	elseif (fcsteqn==2) then
		! 		rf = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow) + kappar(4,1)*log(znow)*log(mnow)) - DELTA
		! 		wf = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rf+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
		! 	else
		! 		wf = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow) + kappaw(4,1)*log(znow)*log(mnow))
		! 		rf = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow) + kappar(4,1)*log(znow)*log(mnow)) - DELTA
		! 	end if
		!
        ! else ! conzflag = 0 or 2
		!
		! 	mp = exp(kappakp(1,1) + kappakp(2,1)*log(mnow) + kappakp(3,1)*log(znow))
		! 	if (fcsteqn==1) then
		! 		wf = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow))
		! 		rf = (znow**(1.0d0/ALPHA))*ALPHA*((wf/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
		! 	elseif (fcsteqn==2) then
		! 		rf = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow)) - DELTA
		! 		wf = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rf+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
		! 	else
		! 		wf = exp(kappaw(1,1) + kappaw(2,1)*log(mnow) + kappaw(3,1)*log(znow))
		! 		rf = exp(kappar(1,1) + kappar(2,1)*log(mnow) + kappar(3,1)*log(znow)) - DELTA
		! 	end if
		!
        ! end if

		call bisectr(agrid,kgrid,bgrid,xgrid,mgrid,znow,mnow,invTa,cfmat0,EVmat0,mumat,ALOW,Yss,idmat,iz,mp, &
			pbij,pxij,time,wf,rf,w0,r0,Lagg0,Ls,Cagg,Hs,Taxagg,Transagg,Tmumat,cdist,hdist,taxdist,transdist,ldist,countmkt)

		wfvec(time) = wf
		w0vec(time) = w0

		Ks = mnow
		Ys = znow*(Ks**ALPHA)*(Lagg0**(1.0d0-ALPHA))
		Gagg = Taxagg-Transagg
		Xagg = Ys-Cagg-Gagg

		! find threshold asset for work decision based on hdist
		!do ix = 1,nx
  !
		!	! ascending order
		!	! NOTE: find minimum ik so that hdist(ik,ix) = 0
  !          do ik = 1,nk
  !              if (hdist(ik,ix) == 0) exit
  !          end do
  !
		!	minasseta(ix) = kgrid(ik)
  !
		!	! descending order
		!	! NOTE: find maximum ik so that hdist(ik,ix) > 0
  !          do ik = nk,1,-1
  !              if (hdist(ik,ix) > 0) exit
  !          end do
  !
		!	! people indexed by kgrid(ik+1) do not work
		!	minassetd(ix) = kgrid(ik+1)
  !
		!end do

		if (mod(time,diagnum)==0) then

	        write(*,"('  time      ', I4, '  ( z, m, w, r) = (', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ')')") &
			time, znow, mnow, w0, r0
			write(*,"('  errors(%)        w : ', F10.5, ', r : ', F10.5, ', inv : ', F10.5)") &
	        log(w0/wf)*100.0d0, log(r0/rf)*100.0d0, log(Xagg/(mp-(1.0d0-DELTA)*mnow))*100.0d0
			write(*,"('  ( w, wf, r, rf) = (', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ') ')") &
			w0, wf, r0, rf
            write(*,"('  bisection ', I4, '  error(%) Lagg-Ls = ', F20.15)") countmkt, abs(log(Lagg0)-log(Ls))*100.0d0
			!do ie = 1,ne
			!	diffma = abs(minasseta(ie)-minassetd(ie))
			!	if (diffma>0.0d0) write(*,"('  WARNING: Multiple thresholds at productivity level ', I4, '  ( ', F8.2, ', ', F8.2, ')')") ie, minasseta(ie), minassetd(ie)
			!end do

			call system_clock(ct2, cr)
		    eptime = dble(ct2-ct1)/cr
		    write(*,"('  Elasped time = ', F10.5)") eptime

		end if

		aggsim(time,1) = Ks
		aggsim(time,2) = w0
		aggsim(time,3) = r0
		aggsim(time,4) = Ys
		aggsim(time,5) = Cagg
		aggsim(time,6) = Xagg
		aggsim(time,7) = Lagg0
		aggsim(time,8) = Hs
		aggsim(time,9) = znow
        aggsim(time,10) = Transagg

        mumat = Tmumat

	end do ! time

	! update bspct
	!bspct = maxval(abs(log(w0vec)-log(wfvec)),1)*2.0d0
    !write(*,"('  updated bspct = ', F8.5)") bspct

    call system_clock(cto2, cro)
    eptime = dble(cto2-cto1)/cro
    write(*,"('  Outerloop elasped time = ', F10.5)") eptime


end subroutine outer


subroutine bisectr(agrid,kgrid,bgrid,xgrid,mgrid,znow,mnow,invTa,cfmat0,EVmat0,mumat,ALOW,Yss,idmat,iz,mp, &
	pbij,pxij,time,wf,rf,w0,r0,Lagg0,Ls,Cagg,Hs,Taxagg,Transagg,Tmumat,cdist,hdist,taxdist,transdist,ldist,countmkt)


	! NOTE: bisectr2: r=rf (or w=wf) is imposed while updating w
	! NOTE: This subroutine returns the equilibrium prices w0 and r0 and quantity Lagg0
	! w0 and r0 are uniquely obtained by bisection (c.f. see the graph of excess demand)
	! Lagg0 is obtained from firm's demand schedule
	! NOTE: 030318 iterationout and maxiterationout are renamed to countmkt and maxcountmkt (these variables are not only for counting # of iterations)
	real(8), intent(in) :: agrid(:), kgrid(:), bgrid(:), xgrid(:), mgrid(:), invTa(:,:), cfmat0(:,:,:,:,:), EVmat0(:,:,:,:), mumat(:,:), znow, mnow, mp, ALOW, Yss
	integer, intent(in) :: idmat(:,:), iz, time !, outermkt, foreind, bsctr, time, maxcountmkt, calcma, bsfixr, linflag
	! real(8), intent(in) :: ALPHA, DELTA, ALOW, B0, HBAR, taul, Yss, Trscale, Transprog1,Transprog2, T0, pbij(:,:), pxij(:,:), tol, tolmkt, bspct
	real(8), intent(in) :: pbij(:,:), pxij(:,:), wf, rf
	real(8), intent(out) :: w0, r0, Lagg0, Ls, Cagg, Hs, Taxagg, Transagg, Tmumat(:,:), cdist(:,:), hdist(:,:), taxdist(:,:), transdist(:,:), ldist(:,:)
	integer, intent(out) :: countmkt
	real(8), parameter :: scale = 10.0d0
	real(8) distancemkt, rL, rH, wL, wH, LaggL, LaggH, g0, LsL, LsH
	integer foreind


	! ! NOTE: will be parameterized
	! ALOW = 0.0d0
	! Yss = 0.0d0

	! NOTE: bisection is done for only r?
	if (time==1) then
		if (calcma==1) print *, "  NOTE: minasset is used to calculate distribution"
		if (outermkt==1) then
			print *, "  NOTE: market-clearing is imposed"
			if (bsctr==1) then
				print *, "  NOTE: bisection is done for interest rate r"
				if (bsfixr==1) print *, "  NOTE: w is fixed at wf while bisecting r"
			else
				print *, "  NOTE: bisection is done for wage w"
				if (bsfixr==1) print *, "  NOTE: r is fixed at rf while bisecting w"
			end if
		end if
	end if

	distancemkt = 1d+4
	countmkt = 0
	foreind = 0 ! ???

	if ((outermkt==0) .or. (foreind==0)) then

		! no market clearing
		countmkt = countmkt + 1
		! subroutine pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,w,r, &
		! 	ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,Ls,Tmumat,cdist,hdist,taxdist,transdist,ldist)
		call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wf,rf, &
			ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,Ls,Tmumat,cdist,hdist,taxdist,transdist,ldist)
		! NOTE: CK assume that excess demand is always equal to zero, hence Lagg (labor demand) = Ls (labor supply).
		Lagg0 = Ls
		w0 = znow*(1.0d0-ALPHA)*(mnow/Lagg0)**ALPHA
		r0 = znow*ALPHA*(mnow/Lagg0)**(ALPHA-1.0d0) - DELTA

	else
	!
	! 	! market clearing by bisection method (using interest rate)
	! 	if (bsctr==1) then
	! 		! NOTE: Initial lower and upper bounds are determined by forecasted interest rate
	! 		rL = (1.0d0-bspct)*rf
	! 		rH = (1.0d0+bspct)*rf
	! 		if (bsfixr==1) then
	! 			wL = wf
	! 			wH = wf
	! 		else
	! 			wL = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rL+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
	! 			wH = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rH+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
	! 		end if
	! 		! Lagg = Firm's FOC implied L
	! 		! r+DELTA = ALPHA*znow*(mnow**(ALPHA-1))*(Lagg**(ALPHA-1))
	! 		LaggL = (( ALPHA*znow*(mnow**(ALPHA-1.0d0)) )/(rL+DELTA) ) ** (1.0d0/(1.0d0-ALPHA))
	! 		LaggH = (( ALPHA*znow*(mnow**(ALPHA-1.0d0)) )/(rH+DELTA) ) ** (1.0d0/(1.0d0-ALPHA))
	!
	! 	else
	!
			wL = (1.0d0-bspct)*wf
			wH = (1.0d0+bspct)*wf
			! if (bsfixr==1) then
			! 	rL = rf
			! 	rH = rf
			! else
				rL = (znow**(1.0d0/ALPHA))*ALPHA*((wL/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
				rH = (znow**(1.0d0/ALPHA))*ALPHA*((wH/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
			! end if
			! Lagg = Firm's FOC implied L
			LaggL = (( (1.0d0-ALPHA)*znow*(mnow**ALPHA) )/wL ) ** (1.0d0/ALPHA)
			LaggH = (( (1.0d0-ALPHA)*znow*(mnow**ALPHA) )/wH ) ** (1.0d0/ALPHA)
	!
	! 	end if
	!
		! Ls = Household implied L
		countmkt = countmkt + 1
		call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wL,rL, &
			ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,LsL,Tmumat,cdist,hdist,taxdist,transdist,ldist)

		countmkt = countmkt + 1
		call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wH,rH, &
			ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,LsH,Tmumat,cdist,hdist,taxdist,transdist,ldist)

		! now countmkt is 2
		! Check if bisection prerequisite satisfies
		do while (( (LsL-LaggL)*(LsH-LaggH)>0.0d0 ) .and. (countmkt<15)) ! try if countmkt = 2, 4, 6, 8, 10, 12, 14

			if (abs(LsL-LaggL) > abs(LsH-LaggH)) then

				print *, "  bisection range is not within market-clearing bound; range extended to right"
	! 			if (bsctr==1) then
	! 				print *, time, rL, rH
	! 			else
					print *, time, wL, wH
	! 			end if
				! write(*,"('  WARNING: Bisection range is not within market-clearing bound; range extended to right ', I4, '  ( ', F8.5, ', ', F8.5, ')')") time, wL, wH
				rL = rH
				wL = wH
				LaggL = LaggH

				countmkt = countmkt + 1
				call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wL,rL, &
					ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,LsL,Tmumat,cdist,hdist,taxdist,transdist,ldist)

	! 			if (bsctr==1) then
	! 				rH = (1.0d0+scale*bspct)*rH
	! 				if (bsfixr==1) then
	! 					wH = wf
	! 				else
	! 					wH = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rH+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
	! 				end if
	! 				LaggH = (( ALPHA*znow*(mnow**(ALPHA-1.0d0)) )/(rH+DELTA) ) ** (1.0d0/(1.0d0-ALPHA))
	! 			else
					wH = (1.0d0+scale*bspct)*wH
					if (bsfixr==1) then
						rH = rf
					else
						rH = (znow**(1.0d0/ALPHA))*ALPHA*((wH/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
					end if
					LaggH = (( (1.0d0-ALPHA)*znow*(mnow**ALPHA) )/wH ) ** (1.0d0/ALPHA)
	! 			end if
	!
				countmkt = countmkt + 1
				call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wH,rH, &
					ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,LsH,Tmumat,cdist,hdist,taxdist,transdist,ldist)

			else

				print *, "  bisection range is not within market-clearing bound; range extended to left"
	! 			if (bsctr==1) then
	! 				print *, time, rL, rH
	! 			else
					print *, time, wL, wH
	! 			end if
				rH = rL
				wH = wL
				LaggH = LaggL

				countmkt = countmkt + 1
				call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wH,rH, &
					ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,LsH,Tmumat,cdist,hdist,taxdist,transdist,ldist)
	!
	! 			if (bsctr==1) then
	! 				rL = (1.0d0-scale*bspct)*rL ! BUG: 180227 rL (not rf) is shifted
	! 				if (bsfixr==1) then
	! 					wL = wf
	! 				else
	! 					wL = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((rL+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
	! 				end if
	! 				LaggL = (( ALPHA*znow*(mnow**(ALPHA-1.0d0)) )/(rL+DELTA) ) ** (1.0d0/(1.0d0-ALPHA))
	! 			else
					wL = (1.0d0-scale*bspct)*wL
					if (bsfixr==1) then
						rL = rf
					else
						rL = (znow**(1.0d0/ALPHA))*ALPHA*((wL/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
					end if
					LaggL = (( (1.0d0-ALPHA)*znow*(mnow**ALPHA) )/wL ) ** (1.0d0/ALPHA)
	! 			end if
	!
				countmkt = countmkt + 1
				call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wL,rL, &
					ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,LsL,Tmumat,cdist,hdist,taxdist,transdist,ldist)

			end if

		end do

	! 	if ((LsL-LaggL)*(LsH-LaggH) > 0.0d0) then
	!
	! 		print *, "  initial range is still not within the market-clearing bound!!"
	! 		print *, "  thus, market-clearing is not imposed"
	! 		countmkt = countmkt + 1
	! 		call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,wf,rf,ALOW,B0,HBAR,taul,Yss,Trscale,Transprog1,Transprog2,T0, &
	! 			pbij,pxij,tol,countmkt,calcma,linflag,Cagg,Hs,Taxagg,Transagg,Ls,Tmumat,cdist,hdist,taxdist,transdist,ldist)
	! 		! CK assumption
	! 		Lagg0 = Ls
	! 		w0 = znow*(1.0d0-ALPHA)*(mnow/Lagg0)**ALPHA
	! 		r0 = znow*ALPHA*(mnow/Lagg0)**(ALPHA-1.0d0) - DELTA
	!
	! 	else
	!
			do while (distancemkt>tolmkt)

	! 			if (bsctr==1) then
	! 				r0 = (rL+rH)/2.0d0
	! 				if (bsfixr==1) then
	! 					w0 = wf
	! 				else
	! 					w0 = (znow**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((r0+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))
	! 				end if
	! 				! Lagg = Firm's FOC implied L
	! 				Lagg0 = (( ALPHA*znow*(mnow**(ALPHA-1.0d0)) )/(r0+DELTA) ) ** (1.0d0/(1.0d0-ALPHA))
	! 			else
					w0 = (wL+wH)/2.0d0
					if (bsfixr==1) then
						r0 = rf
					else
						r0 = (znow**(1.0d0/ALPHA))*ALPHA*((w0/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA
					end if
					! Lagg = Firm's FOC implied L
					Lagg0 = (( (1.0d0-ALPHA)*znow*(mnow**ALPHA) )/w0 ) ** (1.0d0/ALPHA)
	! 			end if
	!
	! 			! Ls = Household implied L
				countmkt = countmkt + 1
				call pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,w0,r0, &
					ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,Ls,Tmumat,cdist,hdist,taxdist,transdist,ldist)
				g0 = (Ls-Lagg0)*(LsL-LaggL) ! f(x0)*f(xl)
	!
	! 			if (bsctr==1) then
	!
	! 				if (g0>0.0d0) then
	! 					rL = r0
	! 				else
	! 					rH = r0
	! 				end if
	! 				distancemkt = abs(log(rH)-log(rL))
	!
	! 			else
	!
					if (g0>0.0d0) then
						wL = w0
					else
						wH = w0
					end if
					distancemkt = abs(log(wH)-log(wL))
	!
	! 			end if
	!
				if (countmkt == maxcountmkt) exit
				! diagnosis
				! write(*,"('  bisection ', I4, '  rH-rL = ', F10.5)") countmkt, distancemkt

			end do ! end of bisection

			r0 = (znow**(1.0d0/ALPHA))*ALPHA*((w0/(1.0d0-ALPHA))**((ALPHA-1.0d0)/ALPHA)) - DELTA

		! end if

	end if


end subroutine bisectr


subroutine pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,w,r, &
	ALOW,Yss,pbij,pxij,Cagg,Hs,Taxagg,Transagg,Ls,Tmumat,cdist,hdist,taxdist,transdist,ldist)
! subroutine pricemap(agrid,kgrid,bgrid,xgrid,mgrid,invTa,cfmat0,EVmat0,mumat,idmat,iz,mp,w,r,ALOW,B0,HBAR,taul,Yss,Trscale,Transprog1,Transprog2,T0, &
! 	pbij,pxij,tol,countmkt,calcma,linflag,Cagg,Hs,Taxagg,Transagg,Ls,Tmumat,cdist,hdist,taxdist,transdist,ldist)


	real(8), intent(in) :: agrid(:), kgrid(:), bgrid(:), xgrid(:), mgrid(:), invTa(:,:), cfmat0(:,:,:,:,:), EVmat0(:,:,:,:), mumat(:,:)
	integer, intent(in) :: idmat(:,:), iz !, countmkt, calcma, linflag
	real(8), intent(in) :: ALOW, Yss, mp, w, r, pbij(:,:), pxij(:,:) !, ALOW, B0, HBAR, taul, Trscale, Transprog1,Transprog2, T0, Yss, pbij(:,:), pxij(:,:), tol
	real(8), intent(out) :: Cagg, Hs, Taxagg, Transagg, Ls, Tmumat(:,:), cdist(:,:), hdist(:,:), taxdist(:,:), transdist(:,:), ldist(:,:)

	real(8), allocatable :: cfmat(:,:,:), EVmat(:,:), EV(:), wkmat(:,:), Wmat(:,:), Nmat(:,:), minasseta(:), minassetd(:)
	integer, allocatable :: jkmat(:,:)
	real(8) betanow, xnow, anow, Trw, Trn, awhigh, anhigh, vwnow, awp, vnnow, anp, ap, wk, atilde, TW, TN, diffma, wm
	integer ia, ik, ie, ib, ix, je, jb, jx, jk, kk, jm


	allocate(EVmat(na,nm),EV(na),jkmat(nk,ne),wkmat(nk,ne),Wmat(na,ne),Nmat(na,ne),minasseta(ne),minassetd(ne))
	if (linflag==2) then
		allocate(cfmat(4,ra+1,1))
	else
		allocate(cfmat(16,ra+1,rm+1))
	end if

    if (Trscale == 0.0d0) then
        ! for a below atilde, not-working is not feasible (consumption negative)
        atilde = (alow)/(1.0d0+r)
    else
        call transfers(0.0d0,0.0d0,r,w,0.0d0,Yss,Trscale,Transprog,T0,Trn)
        atilde = (alow - Trn)/(1.0d0+r)
    end if


	if (calcma==1) then
		! ****************************************
		! threshold asset for work decision
		minasseta = 0.0d0 ! ascending order
		minassetd = 0.0d0 ! descending order

		!$omp parallel do default(shared) private(ie,ib,ix,betanow,xnow,ia,anow,awp,anp,TW,TN,EVmat,cfmat,jm,wm,EV)
		do ie = 1,ne

			ib = idmat(ie,1)
			ix = idmat(ie,2)
			betanow = bgrid(ib)
			xnow = xgrid(ix)

			! now we know ie
			EVmat = EVmat0(:,:,ie,iz)
			if (linflag==0) then
				cfmat = cfmat0(:,:,:,ie,iz)
			! elseif (linflag==1) then
			! 	EVmat = EVmat0(:,:,ie,iz)
			elseif (linflag==2) then
				call findiw(mp,mgrid,jm,wm)
				EV = wm*EVmat0(:,jm,ie,iz) + (1.0d0-wm)*EVmat0(:,jm+1,ie,iz)
				cfmat(:,:,1) = spfit(invTa,EV,ra,agrid)
			end if

			do ia = 1,na

				anow = agrid(ia)
				call calcwnval(xnow,anow,r,w,mp,ALOW,Yss,atilde,betanow,agrid,mgrid,cfmat,EVmat,awp,anp,TW,TN,0)
				Wmat(ia,ie) = TW
				Nmat(ia,ie) = TN

			end do

		end do
		!$omp end parallel do

		!$omp parallel do default(shared) private(ie,ib,ix,betanow,xnow,EVmat,cfmat,jm,wm,EV)
		do ie = 1,ne

			ib = idmat(ie,1)
			ix = idmat(ie,2)
			betanow = bgrid(ib)
			xnow = xgrid(ix)

			! now we know ie
			EVmat = EVmat0(:,:,ie,iz)
			if (linflag==0) then
				cfmat = cfmat0(:,:,:,ie,iz)
			elseif (linflag==2) then
				call findiw(mp,mgrid,jm,wm)
				EV = wm*EVmat0(:,jm,ie,iz) + (1.0d0-wm)*EVmat0(:,jm+1,ie,iz)
				cfmat(:,:,1) = spfit(invTa,EV,ra,agrid)
			end if

			minasseta(ie) = bisectma(agrid,mgrid,cfmat,EVmat,Wmat(:,ie),Nmat(:,ie),xnow,r,w,mp,Yss,atilde,HBAR,Trscale,Transprog,T0,betanow,B0,taul,ALOW,tol,1,linflag)
			minassetd(ie) = bisectma(agrid,mgrid,cfmat,EVmat,Wmat(:,ie),Nmat(:,ie),xnow,r,w,mp,Yss,atilde,HBAR,Trscale,Transprog,T0,betanow,B0,taul,ALOW,tol,0,linflag)

		end do
		!$omp end parallel do

		if (sum(abs(minasseta-minassetd),1)>tol) then
            ! write(*,"('  multiple thresholds at bisection ', I4)") countmkt
			! write(*,"('  multiple thresholds at bisection ', I4)") countmkt
		 	do ie = 1,ne
		 		diffma = abs(minasseta(ie)-minassetd(ie))
		 		if (diffma>0) print *, ie, minasseta(ie), minassetd(ie)
            end do
        end if

	end if


	hdist = 0.0d0
	cdist = 0.0d0
	taxdist = 0.0d0
	transdist = 0.0d0
	ldist = 0.0d0

	!$omp parallel do default(shared) private(ie,ib,ix,betanow,xnow,ik,anow,Trw,Trn,awhigh,anhigh,vwnow,awp,vnnow,anp,ap,jk,wk,kk,EVmat,cfmat,jm,wm,EV)
	do ik = 1,nk

		anow = kgrid(ik)

		do ie = 1,ne

			ib = idmat(ie,1)
			ix = idmat(ie,2)
			betanow = bgrid(ib)
			xnow = xgrid(ix)

			! now we know ie
			EVmat = EVmat0(:,:,ie,iz)
			if (linflag==0) then
				cfmat = cfmat0(:,:,:,ie,iz)
			elseif (linflag==2) then
				call findiw(mp,mgrid,jm,wm)
				EV = wm*EVmat0(:,jm,ie,iz) + (1.0d0-wm)*EVmat0(:,jm+1,ie,iz)
				cfmat(:,:,1) = spfit(invTa,EV,ra,agrid)
            end if

			! NOTE: about 80% of hshlds are active
            ! 20181121 bug: below if option makes multiple threshold warning active for irrelevant reason
			if (mumat(ik,ie)>0.0d0) then

				if (calcma==1) then
					! NOTE: if we know minasset beforehand, # of calling gss is nk*ne + 2*na*ne + (# in bisection) < 2*nk*ne
					if (anow >= minassetd(ie)) then
						! Not work
						call transfers(xnow,anow,r,w,0.0d0,Yss,Trscale,Transprog,T0,Trn)
						anhigh = (1.0d0+r)*anow + Trn
						call gss(ALOW,anhigh,anp,tol**2,0,linflag,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Trn,betanow,B0,Hbar,taul)
						ap = anp
						hdist(ik,ie) = 0.0d0
						cdist(ik,ie) = anhigh-ap
	                    taxdist(ik,ie) = 0.0d0
	                    transdist(ik,ie) = Trn
					else
						! Work
						call transfers(xnow,anow,r,w,HBAR,Yss,Trscale,Transprog,T0,Trw)
						awhigh = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Trw
						call gss(ALOW,awhigh,awp,tol**2,1,linflag,anow,xnow,agrid,mgrid,cfmat,EVmat,mp,w,r,Trw,betanow,B0,Hbar,taul)
						ap = awp
						hdist(ik,ie) = HBAR
						cdist(ik,ie) = awhigh-ap
	                    taxdist(ik,ie) = taul*w*xnow*HBAR
	                    transdist(ik,ie) = Trw
					end if

				else
					! NOTE: We calculate vw and vn and then compare these values... may cause approximation errors? (i.e. nonmonotone hdist)
					! # of calling gss is 2*nk*ne
					call calcwnval(xnow,anow,r,w,mp,ALOW,Yss,atilde,betanow,agrid,mgrid,cfmat,EVmat,awp,anp,vwnow,vnnow,0)

					if (vwnow>vnnow) then
						call transfers(xnow,anow,r,w,HBAR,Yss,Trscale,Transprog,T0,Trw)
						awhigh = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Trw
						ap = awp
						hdist(ik,ie) = HBAR
						cdist(ik,ie) = awhigh-ap
	                    taxdist(ik,ie) = taul*w*xnow*HBAR
	                    transdist(ik,ie) = Trw
					else
						call transfers(xnow,anow,r,w,0.0d0,Yss,Trscale,Transprog,T0,Trn)
						anhigh = (1.0d0+r)*anow + Trn
						ap = anp
						hdist(ik,ie) = 0.0d0
						cdist(ik,ie) = anhigh-ap
	                    taxdist(ik,ie) = 0.0d0
	                    transdist(ik,ie) = Trn
	                end if

				end if

                ldist(ik,ie) = xnow*hdist(ik,ie)

				! finding indices and weights
				! NOTE: finding indices and weights on grid points of k in the next periods
				! Ref. Young (2010?) and Rios-Rull (1999?)
				call findjw(ap,kgrid,jk,wk)
				jkmat(ik,ie) = min(max(1,jk),nk-1)
				wkmat(ik,ie) = wk

			end if

		end do

	end do
	!$omp end parallel do

	! transition
	Tmumat = 0.0d0
	do ie = 1,ne

		ib = idmat(ie,1)
		ix = idmat(ie,2)

		do ik = 1,nk

			if (mumat(ik,ie) > 0.0d0) then

				jk = jkmat(ik,ie)
				wk = wkmat(ik,ie)

				do je = 1,ne

					jb = idmat(je,1)
					jx = idmat(je,2)
					Tmumat(jk,je) = Tmumat(jk,je) + pbij(ib,jb)*pxij(ix,jx)*wk*mumat(ik,ie)
					Tmumat(jk+1,je) = Tmumat(jk+1,je) + pbij(ib,jb)*pxij(ix,jx)*(1.0d0-wk)*mumat(ik,ie)

				end do

			end if

		end do

	end do

	Cagg = 0.0d0
	Hs = 0.0d0
	Taxagg = 0.0d0
	Transagg = 0.0d0
	Ls = 0.0d0

	do ie = 1,ne

		Cagg = Cagg + sum(mumat(:,ie)*cdist(:,ie),1)
		Hs   = Hs   + sum(mumat(:,ie)*hdist(:,ie),1)
		Taxagg = Taxagg + sum(mumat(:,ie)*taxdist(:,ie),1)
		Transagg = Transagg + sum(mumat(:,ie)*transdist(:,ie),1)
		Ls   = Ls   + sum(mumat(:,ie)*ldist(:,ie),1)

	end do


end subroutine pricemap


function bisectma(agrid,mgrid,cfmat,EVmat,Wvec,Nvec,xnow,r,w,mp,Yagg,atilde,HBAR,Trscale,Transprog,T0,betanow,B0,taul,ALOW,tol,ascend,linflag) result(ma0)


	real(8), intent(in) :: agrid(:), mgrid(:), cfmat(:,:,:), EVmat(:,:), Wvec(:), Nvec(:), xnow, r, w, mp, Yagg, atilde, &
		HBAR, Trscale, Transprog, T0, betanow, B0, taul, ALOW, tol
	integer, intent(in) :: ascend, linflag
	real(8) ma0, malow, mahigh, WNdiffL, WNdiffH, WNdiff0, awp, anp, TW, TN, distance
	integer ia, na, iteration


	na = size(agrid,1)

	if (ascend) then
		! ascending order
		! NOTE: find minimum ia so that W(ia)<N(ia)
		do ia = 1,na
		    if ((Wvec(ia)-Nvec(ia)) < 0.0d0) exit
		end do

	else
		! descending order
		! NOTE: find minimum ia so that W(ia)>N(ia)
		do ia = na,1,-1
		    if ((Wvec(ia)-Nvec(ia)) > 0.0d0) exit
		end do
		ia = ia + 1

	end if

	if (ia == 1) then

	    ma0 = agrid(1)

	elseif (ia > na) then

	    ma0 = agrid(na)

	else

	    malow = agrid(ia-1)
	    WNdiffL = Wvec(ia-1) - Nvec(ia-1)
	    mahigh = agrid(ia)
	    WNdiffH = Wvec(ia) - Nvec(ia)

	    if (WNdiffL*WNdiffH > 0) then

		    print *, "  Bisection prerequisite violated!"
	        write(*,"('  ia-1 ', I4, '   ia ', I4)") ia-1, ia

		else

			! bisection
			distance = 1d+4
	        iteration = 0

	        do while (distance > tol*10.0d0)

	            ma0 = (malow+mahigh)/2.0d0

				call calcwnval(xnow,ma0,r,w,mp,ALOW,Yagg,atilde,betanow,agrid,mgrid,cfmat,EVmat,awp,anp,TW,TN,0)

	            WNdiff0 = TW-TN

	            ! f(x0)*f(xl)
	            if (WNdiff0*WNdiffL > 0) then
	                malow = ma0
	            else
	                mahigh = ma0
	            end if

	            distance = mahigh-malow
	            iteration = iteration + 1
	            ! diagnosis
	            !write(*,"('  bisection ', I4, '  ahigh-alow = ', F10.5)") iteration, distance

	        end do ! while loop for bisection

	    end if

	end if


end function bisectma


subroutine calcforecast(drop,izvec,zgrid,ymat,DELTA,kappakpnew,kappawnew,kapparnew,kappalnew,rsq,conzflag)

	integer, intent(in) :: drop, izvec(:), conzflag
	real(8), intent(in) :: ymat(:,:), DELTA, zgrid(:)
	real(8), intent(out) :: kappakpnew(:,:), kappawnew(:,:), kapparnew(:,:), kappalnew(:,:), rsq(:,:)
	integer iz, nz, simTT, num
    real(8), allocatable :: mvec(:), wvec(:), rvec(:), lvec(:), XX(:,:), yy(:,:), ee(:,:), yyhat(:,:), zvec(:)
    real(8) AA(2,2), bb(2,4), det, invAA(2,2), kappa(2,4), temp(1), yymean(4), AA0(3,3), bb0(3,4), invAA0(3,3), kappa0(3,4)
    real(8) AA3(4,4), bb3(4,4), invAA3(4,4), kappa3(4,4)   ! for option: conzflag == 3


	! 23 May 2018: kappakp shape has been altered.
    !nz = size(kappakpnew,1)
    nz = size(kappakpnew,2)
	simTT = size(ymat,1)
	allocate(mvec(simTT),wvec(simTT),rvec(simTT),lvec(simTT),zvec(simTT))
	mvec = ymat(:,1)
	wvec = ymat(:,2)
    ! rvec is the rental rate
	rvec = ymat(:,3) + DELTA
    lvec = ymat(:,7)
	zvec = ymat(:,9)

	if (conzflag==1) then

        do iz = 1,nz

		    num = count(izvec(drop+1:simTT-1)==iz)

		    allocate(XX(num,2),yy(num,4),ee(num,4),yyhat(num,4))
		    XX(:,1) = 1.0d0
		    XX(:,2) = log(pack(mvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))
		    yy(:,1) = log(pack(mvec(drop+2:simTT),  izvec(drop+1:simTT-1)==iz))
		    yy(:,2) = log(pack(wvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))
		    yy(:,3) = log(pack(rvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))
            yy(:,4) = log(pack(lvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))

            yymean = sum(yy,1)/num

		    AA = matmul(transpose(XX),XX)
		    bb = matmul(transpose(XX),yy)
		    det = AA(2,2)*AA(1,1) - AA(1,2)*AA(2,1)
		    invAA(1,:) = (/AA(2,2), -AA(1,2)/)
		    invAA(2,:) = (/-AA(2,1), AA(1,1)/)
		    invAA =  invAA/det
		    kappa = matmul(invAA,bb)

			kappakpnew(1,iz) = kappa(1,1)
			kappakpnew(2,iz) = kappa(2,1)
		    kappawnew(1,iz) = kappa(1,2)
			kappawnew(2,iz) = kappa(2,2)
		    kapparnew(1,iz) = kappa(1,3)
			kapparnew(2,iz) = kappa(2,3)
            kappalnew(1,iz) = kappa(1,4)
			kappalnew(2,iz) = kappa(2,4)

		    yyhat(:,1) = kappa(1,1) + kappa(2,1)*XX(:,2)
		    yyhat(:,2) = kappa(1,2) + kappa(2,2)*XX(:,2)
		    yyhat(:,3) = kappa(1,3) + kappa(2,3)*XX(:,2)
            yyhat(:,4) = kappa(1,4) + kappa(2,4)*XX(:,2)

			ee = yy - yyhat
		    temp = matmul(reshape(ee(:,1),(/1,num/)),ee(:,1))
		    temp = temp/matmul(reshape(yy(:,1)-yymean(1),(/1,num/)),yy(:,1)-yymean(1))
		    rsq(1,iz) = 1.0d0-temp(1)
		    temp = matmul(reshape(ee(:,2),(/1,num/)),ee(:,2))
		    temp = temp/matmul(reshape(yy(:,2)-yymean(2),(/1,num/)),yy(:,2)-yymean(2))
		    rsq(2,iz) = 1.0d0-temp(1)
		    temp = matmul(reshape(ee(:,3),(/1,num/)),ee(:,3))
		    temp = temp/matmul(reshape(yy(:,3)-yymean(3),(/1,num/)),yy(:,3)-yymean(3))
		    rsq(3,iz) = 1.0d0-temp(1)
            temp = matmul(reshape(ee(:,4),(/1,num/)),ee(:,4))
		    temp = temp/matmul(reshape(yy(:,4)-yymean(4),(/1,num/)),yy(:,4)-yymean(4))
		    !temp = temp/matmul(reshape(yyhat(:,4),(/1,num/)),yyhat(:,4))
            ! BUG found on 8 March 2018
		    rsq(4,iz) = 1.0d0-temp(1)

		    deallocate(XX,yy,ee,yyhat)

        end do

    elseif (conzflag == 3) then

        num = simTT - drop - 1

        allocate(XX(num,4),yy(num,4),ee(num,4),yyhat(num,4))
		XX(:,1) = 1.0d0
		XX(:,2) = log(mvec(drop+1:simTT-1))
        XX(:,3) = log(zvec(drop+1:simTT-1))
        XX(:,4) = log(zvec(drop+1:simTT-1))*log(mvec(drop+1:simTT-1))
		yy(:,1) = log(mvec(drop+2:simTT))
		yy(:,2) = log(wvec(drop+1:simTT-1))
		yy(:,3) = log(rvec(drop+1:simTT-1))
        yy(:,4) = log(lvec(drop+1:simTT-1))


        yymean = sum(yy,1)/num

		AA3 = matmul(transpose(XX),XX)
		bb3 = matmul(transpose(XX),yy)
		! NOTE: using mkl-lapack
        invAA3 = inv(AA3)
		kappa3 = matmul(invAA3,bb3)

		kappakpnew(1,1) = kappa3(1,1)
        kappakpnew(2,1) = kappa3(2,1)
        kappakpnew(3,1) = kappa3(3,1)
        kappakpnew(4,1) = kappa3(4,1)
        kappawnew(1,1) = kappa3(1,2)
        kappawnew(2,1) = kappa3(2,2)
        kappawnew(3,1) = kappa3(3,2)
        kappawnew(4,1) = kappa3(4,2)
		kapparnew(1,1) = kappa3(1,3)
        kapparnew(2,1) = kappa3(2,3)
        kapparnew(3,1) = kappa3(3,3)
        kapparnew(4,1) = kappa3(4,3)
        kappalnew(1,1) = kappa3(1,4)
        kappalnew(2,1) = kappa3(2,4)
        kappalnew(3,1) = kappa3(3,4)
        kappalnew(4,1) = kappa3(4,4)

		yyhat(:,1) = kappa3(1,1) + kappa3(2,1)*XX(:,2) + kappa3(3,1)*XX(:,3) + kappa3(4,1)*XX(:,4)
		yyhat(:,2) = kappa3(1,2) + kappa3(2,2)*XX(:,2) + kappa3(3,2)*XX(:,3) + kappa3(4,2)*XX(:,4)
		yyhat(:,3) = kappa3(1,3) + kappa3(2,3)*XX(:,2) + kappa3(3,3)*XX(:,3) + kappa3(4,3)*XX(:,4)
        yyhat(:,4) = kappa3(1,4) + kappa3(2,4)*XX(:,2) + kappa3(3,4)*XX(:,3) + kappa3(4,4)*XX(:,4)
		ee = yy - yyhat
		temp = matmul(reshape(ee(:,1),(/1,num/)),ee(:,1))
		temp = temp/matmul(reshape(yy(:,1)-yymean(1),(/1,num/)),yy(:,1)-yymean(1))
		rsq(1,1) = 1.0d0-temp(1)
		temp = matmul(reshape(ee(:,2),(/1,num/)),ee(:,2))
		temp = temp/matmul(reshape(yy(:,2)-yymean(2),(/1,num/)),yy(:,2)-yymean(2))
		rsq(2,1) = 1.0d0-temp(1)
		temp = matmul(reshape(ee(:,3),(/1,num/)),ee(:,3))
		temp = temp/matmul(reshape(yy(:,3)-yymean(3),(/1,num/)),yy(:,3)-yymean(3))
		rsq(3,1) = 1.0d0-temp(1)
        temp = matmul(reshape(ee(:,4),(/1,num/)),ee(:,4))
		temp = temp/matmul(reshape(yy(:,4)-yymean(4),(/1,num/)),yy(:,4)-yymean(4))
		!temp = temp/matmul(reshape(yyhat(:,4),(/1,num/)),yyhat(:,4))
        ! BUG found on 8 March 2018
		rsq(4,1) = 1.0d0-temp(1)

		deallocate(XX,yy,ee,yyhat)

	else

        num = simTT - drop - 1

        allocate(XX(num,3),yy(num,4),ee(num,4),yyhat(num,4))
		XX(:,1) = 1.0d0
		XX(:,2) = log(mvec(drop+1:simTT-1))
        XX(:,3) = log(zvec(drop+1:simTT-1))
		yy(:,1) = log(mvec(drop+2:simTT))
		yy(:,2) = log(wvec(drop+1:simTT-1))
		yy(:,3) = log(rvec(drop+1:simTT-1))
        yy(:,4) = log(lvec(drop+1:simTT-1))

        yymean = sum(yy,1)/num

		AA0 = matmul(transpose(XX),XX)
		bb0 = matmul(transpose(XX),yy)
		! NOTE: using matrix inverse formula
		! det = AA0(3,3)*AA0(2,2)*AA(1,1) + AA0(1,2)*AA0(2,3)*AA0(3,1) + AA0(1,3)*AA0(2,1)*AA0(3,2) - AA0(1,3)*AA0(2,2)*AA0(3,1) - AA0(1,1)*AA0(2,3)*AA0(3,2) - AA0(1,2)*AA0(2,1)*AA0(3,3)
		! invAA0(1,1) = AA0(2,2)*AA0(3,3) - AA0(2,3)*AA0(3,2)
		! invAA0(1,2) = AA0(1,3)*AA0(3,2) - AA0(1,2)*AA0(3,3)
		! invAA0(1,3) = AA0(1,2)*AA0(2,3) - AA0(1,3)*AA0(2,2)
		! invAA0(2,1) = AA0(2,3)*AA0(3,1) - AA0(2,1)*AA0(3,3)
		! invAA0(2,2) = AA0(1,1)*AA0(3,3) - AA0(1,3)*AA0(3,1)
		! invAA0(2,3) = AA0(1,3)*AA0(2,1) - AA0(1,1)*AA0(2,3)
		! invAA0(3,1) = AA0(2,1)*AA0(3,2) - AA0(2,2)*AA0(3,1)
		! invAA0(3,2) = AA0(1,2)*AA0(3,1) - AA0(1,1)*AA0(3,2)
		! invAA0(3,3) = AA0(1,1)*AA0(2,2) - AA0(1,2)*AA0(2,1)
		!
		! invAA0 = invAA0/det
		! NOTE: using mkl-lapack
        invAA0 = inv(AA0)
		kappa0 = matmul(invAA0,bb0)

		kappakpnew(1,1) = kappa0(1,1)
        kappakpnew(2,1) = kappa0(2,1)
        kappakpnew(3,1) = kappa0(3,1)
        kappawnew(1,1) = kappa0(1,2)
        kappawnew(2,1) = kappa0(2,2)
        kappawnew(3,1) = kappa0(3,2)
		kapparnew(1,1) = kappa0(1,3)
        kapparnew(2,1) = kappa0(2,3)
        kapparnew(3,1) = kappa0(3,3)
        kappalnew(1,1) = kappa0(1,4)
        kappalnew(2,1) = kappa0(2,4)
        kappalnew(3,1) = kappa0(3,4)

		yyhat(:,1) = kappa0(1,1) + kappa0(2,1)*XX(:,2) + kappa0(3,1)*XX(:,3)
		yyhat(:,2) = kappa0(1,2) + kappa0(2,2)*XX(:,2) + kappa0(3,2)*XX(:,3)
		yyhat(:,3) = kappa0(1,3) + kappa0(2,3)*XX(:,2) + kappa0(3,3)*XX(:,3)
        yyhat(:,4) = kappa0(1,4) + kappa0(2,4)*XX(:,2) + kappa0(3,4)*XX(:,3)
		ee = yy - yyhat
		temp = matmul(reshape(ee(:,1),(/1,num/)),ee(:,1))
		temp = temp/matmul(reshape(yy(:,1)-yymean(1),(/1,num/)),yy(:,1)-yymean(1))
		rsq(1,1) = 1.0d0-temp(1)
		temp = matmul(reshape(ee(:,2),(/1,num/)),ee(:,2))
		temp = temp/matmul(reshape(yy(:,2)-yymean(2),(/1,num/)),yy(:,2)-yymean(2))
		rsq(2,1) = 1.0d0-temp(1)
		temp = matmul(reshape(ee(:,3),(/1,num/)),ee(:,3))
		temp = temp/matmul(reshape(yy(:,3)-yymean(3),(/1,num/)),yy(:,3)-yymean(3))
		rsq(3,1) = 1.0d0-temp(1)
        temp = matmul(reshape(ee(:,4),(/1,num/)),ee(:,4))
		temp = temp/matmul(reshape(yy(:,4)-yymean(4),(/1,num/)),yy(:,4)-yymean(4))
		!temp = temp/matmul(reshape(yyhat(:,4),(/1,num/)),yyhat(:,4))
        ! BUG found on 8 March 2018
		rsq(4,1) = 1.0d0-temp(1)

		deallocate(XX,yy,ee,yyhat)

	end if


end subroutine calcforecast


end module mod_outer

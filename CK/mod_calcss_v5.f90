module mod_calcss
! v2 incorporates Howard's improvement algorithm 4 May 2018
! v3 incorporates simulation
! v4 considers a different approach to find market clearing
! v5 uses transfers by wealth

use mod_randnumber
use mod_utils
use mod_spline
use mkl_spblas
use mod_policies
use mod_parameters


implicit none


contains


! subroutine calcss(knotsk,knotsb,invTk,znow,lnow,Ge,Gy,Gd,Pe,Py,Pd,idmat,muxz,vmat0,gmat0,mu0,evec,mpmat)
subroutine calcss(agrid,kgrid,bgrid,xgrid,pbij,pxij,xss,ALOW,idmat,calcdist,calcdistsim, &
	    w,r,zagg,Lagg,Kagg,Tragg,mumat,ndist,Wmat,Nmat,minasset,Ls,Ks,EPratio,TrW5,EmpW5,stdlogwage,stdlogearn,rhoxest,error)

    ! use mod_parameters
	integer, parameter :: savessresult = 1
	! integer, parameter :: calcdist = 1
	! integer, parameter :: calcdistsim = 1

	real(8), intent(out) :: agrid(:), kgrid(:), bgrid(:), pbij(:,:), xgrid(:), pxij(:,:), xss(:), ALOW
	integer, intent(in)  :: idmat(:,:), calcdist, calcdistsim
	real(8), intent(in)  :: zagg
	real(8), intent(out) :: w, r, Lagg, Kagg, Tragg, mumat(:,:), ndist(:,:), Wmat(:,:), Nmat(:,:), minasset(:)
    real(8), intent(out) :: Ls, Ks, EPratio, TrW5(:), EmpW5(:), stdlogwage, stdlogearn, rhoxest
    integer, intent(out) :: error

	integer howard, howardcnt
	real(8) Vmat(na,ne) !, Wmat(na,ne), Nmat(na,ne)
	real(8) apmat(na,ne), awpmat(na,ne), anpmat(na,ne), lmat(na,ne)
    real(8) apmatnew(na,ne), awpmatnew(na,ne), anpmatnew(na,ne), lmatnew(na,ne)
	real(8) TVmat(na,ne), TWmat(na,ne), TNmat(na,ne)
	integer jkmat(nk,ne)
	real(8) EV(na), Tmumat(nk,ne), wkmat(nk,ne), muk(nk) !, minasset(ne)
    real(8) awpmat1(nk*10,ne), anpmat1(nk*10,ne), kgrid1(nk*10)
	real(8) AA(nk*ne,nk*ne), muvec(nk*ne)

    integer ib, jb, ie, je, ix, jx, ia, ik, jk, kk, iteration, iterationout, i
    real(8) Yagg, TrY, khigh, khighnew, Kratio, betanow, xnow, anow, ap, atilde, kwhigh, knhigh, awp, anp, TW, TN, earns, tax, Tr, Trw, Trn, wk
    real(8) distance, distanceout, damp, wnmax, wnmin, distancevec(maxiter), distancedamp
    real(8) rnew, wnew
    ! integer, allocatable :: jkmat(:,:)
	! real(8), allocatable :: EV(:), Tmumat(:,:), wkmat(:,:), muk(:), minasset(:)

    real(8) musum, Ys, Govt, KY, TaxY, Taxagg, Giniwage, Giniearns, Giniwealth, distpolmax(3), distpolmean(3), TrInc5(5), EPsim
    real(8) ShareW51, ShareW52, ShareW53, ShareW54, ShareW55, EmpW51, EmpW52, EmpW53, EmpW54, EmpW55, TrInc51, TrInc52, TrInc53, TrInc54, TrInc55, TrW51, TrW52, TrW53, TrW54, TrW55, objval

	! ! for SimulSS_v4b
    ! integer, parameter :: NS = 100000 ! number of samples
    ! integer, parameter :: T = 100 ! number of periods (quarters) should be 4x
    ! real(8) distind(NS), temp(NS), asim(NS,T+1)
    ! integer xsimind(NS,T), xsimindnow(NS), xsimindnext(NS), xsimtemp(1), xindnow(NS)
    ! integer, allocatable :: xindnext(:,:)
    ! integer asimindtemp, maxsimind
    ! ! real(8) minasset(ne), awpmat1(nk*10,ne), anpmat1(nk*10,ne), k1grid(nk*10)
    ! real(8) worksim(NS,T), consim(NS,T), earningsim(NS,T), wagesim(NS,T), wage0sim(NS,T), incomesim(NS,T), postincomesim(NS,T), Trsim(NS,T)
    ! real(8) Et(T), EPt(T)
    ! ! real(8) atilde, anow, xnow, Tr, EPsim, stdlogwage, stdlogearn
	! ! real(8) EPsim, stdlogwage, stdlogearn
    ! real(8) anwagesim(NS,int(T/4)), anearnsim(NS,int(T/4)), fulltime(NS,int(T/4)), logwage(NS,int(T/4)), logearn(NS,int(T/4))
    ! real(8) lnwagetemp(int(NS*T/4)), lnwagelagtemp(int(NS*T/4))
    ! real(8), allocatable :: mumeasure(:), xindcount(:), logwage_active(:), logearn_active(:), lnwagesim(:), lnwagelagsim(:), X(:,:)
    ! real(8) XX(2,2), Xy(2), estimate(2) !, rhoxest
    ! integer Nactive, count_active, time

    ! for elapsed time
	integer cto1, cto2, cro
    real(8) eptimeo
    integer ct1, ct2, cr
    real(8) eptime
    ! for vsl
    integer errcode, brng, method, seed
    type(vsl_stream_state) :: stream
	! for dgetrf, dgetri
    integer ra, INFO
    integer, allocatable :: IPIV(:)
    real(8), allocatable :: invTa(:,:), cfa(:,:), WORK(:)
	! for mkl_spblas
	integer csrColInd(nk*ne*ne*2)
    ! integer, allocatable :: csrColInd(:)
    integer csrRowPtr(nk*ne+1)
	real(8) csrVal(nk*ne*ne*2)
    ! real(8), allocatable :: csrVal(:)
    ! Matrix descriptor
    TYPE(MATRIX_DESCR) descrA     ! Sparse matrix descriptor
    ! CSR matrix representation
    TYPE(SPARSE_MATRIX_T) csrA    ! Structure with sparse matrix
    integer index !, i, info


	call system_clock(cto1, cro)

	ra = na-2
    allocate(invTa(ra,ra),cfa(4,ra+1),IPIV(ra),WORK(ra))


	! NOTE: Computed here for calibration???
	! grids for preference shocks a la Krusell and Smith (1998): to be done
	pbij = 1.0d0
    bgrid = BETA
	! grids for labor productivity
	if (tauflag==1) then
		call tauchen(nx,0.0d0,RHOX,SDINOVX,mx,xgrid,pxij)
	elseif (tauflag==0) then
		call rouwenhorst(RHOX,SDINOVX,nx,xgrid,pxij)
	end if
    pxij = max(pxij,0.0d0)
    xgrid = exp(xgrid)
    call markovss(pxij,nx,tol,xss)
    ! write(*,"('  xss =', F10.5)") xss


    ! ********************************************************************************
    !  OUTER LOOP: find equilibrium prices
    ! ********************************************************************************
    ! initial guess for L (not important)
    Lagg = 0.30d0
    ! 8 May 2018: initial guess for r
    r = 0.01d0
    Kagg = Lagg*( ((r + DELTA)/ALPHA/zagg )**(1.0d0/(ALPHA-1.0d0)) )
    w = (1.0d0-ALPHA)*zagg*(Kagg/Lagg)**ALPHA

    ! initial guess for the value function
	! NOTE: The value function converged in the previous inner loop is used as an initial guess in the current inner loop
	! We need an initial guess for the value function for the first time
    Vmat = 0.0d0

	! initial distribution
    ! NOTE: 22 Jan 2018 The initial distribution below is much more stable
	mumat = 1.0d0/dble(nk*ne)
 	! 26 Feb 2018: needed for values in the other cells of mumat, otherwise random values are set
	! mumat = 0.0d0
	! do ie = 1,ne
    !     mumat(11:20,ie) = 0.1d0/dble(ne)
    ! end do

    distanceout = 1d+4
    iterationout = 1
    damp = dampss ! NOTE: needed?


    error = 0

    do while (distanceout > tolout)
	! do while (iterationout<2)

        ! NOTE: given r and Lagg, Kagg and w consistent with firm FOC
        ! NOTE: Kagg,  and Yagg are capital, labor and output implied by the demand (firm) side
        ! Kagg = *( ((r + DELTA)/ALPHA/zagg )**(1.0d0/(ALPHA-1.0d0)) )
        ! w = (zagg**(1.0d0/(1.0d0-ALPHA)))*(1.0d0-ALPHA)*(((r+DELTA)/ALPHA)**(ALPHA/(ALPHA-1.0d0)))

        ! for borrowing constraint (not used)
        ! Yagg = zagg*(Kagg**(ALPHA))*(Lagg**(1.d0-ALPHA))
        ! ALOW = -phi*Yagg
        ALOW = 0.0d0

		! for a below atilde, not-working is not feasible (consumption negative)
        if (Trscale == 0.0d0) then
            atilde = (alow)/(1.0d0+r)
        else
			! maximum transfer to those with no wealth
            call transfers(0.0d0,0.0d0,r,w,0.0d0,Yagg,Trscale,Transprog,T0,Trn)
            atilde = (alow - Trn)/(1.0d0+r)
        end if

		! khigh = 250.0d0 ! for nx=5
		khigh = 700.0d0 ! for nx=7 or 11
		! khigh = 2000.0d0 ! for nx=17

        ! adaptive maximum of grids
        ! if (iterationout == 1) then
        !     khigh = Kagg*200.0d0
        ! else
		!
		! 	! NOTE: find maximum ik so that muk(ik)>0
        !     do ik = 1,nk
        !         if (muk(nk-ik+1)>0.0d0) exit
        !     end do
		!
        !     khighnew = 0.7d0*khigh + (1-0.3d0)*kgrid(nk-ik+1)*1.3d0
		!
        !     if  ((abs(log(khighnew)-log(khigh)) > 0.05d0) .and. (iterationout < 20)) then
        !         khigh = khighnew
        !     end if
		!
		! end if
		!
        ! Kratio = khigh/Kagg

		write(*,*) " "
		print *,'===================================================================='
		print *,' iteration ', iterationout
		write(*,"('  r =', F10.5, '    w =', F10.5, '   blimit =', F10.5)") r, w, alow
        ! write(*,"('  khigh over K =', F10.2)") Kratio
		write(*,"('  khigh =', F10.2)") khigh

		! grids for assets
		if (congrid==1) then
			! more grid points concentrated toward the left
			agrid = logspace(0.0d0,deggrid,na) - 10.0d0**(0.0d0)
			agrid = agrid/(10.0d0**deggrid - 10.0d0**0.0d0)
			agrid = agrid*(khigh-alow) + alow
			! kgrid = logspace(0.0d0,deggrid,nk) - 10.0d0**(0.0d0)
			! kgrid = kgrid/(10.0d0**deggrid - 10.0d0**0.0d0)
			! kgrid = kgrid*(khigh-alow) + alow
            kgrid = logspace(log(ALOW - 1.0d0*ALOW + 1.0d0)/log(10.0d0), log(khigh - 1.0d0*ALOW + 1.0d0)/log(10.0d0), nk)
	        kgrid = kgrid + ALOW - 1.0d0

		else
	        agrid = logspace(log(ALOW - 1.0d0*ALOW + 1.0d0)/log(10.0d0), log(khigh - 1.0d0*ALOW + 1.0d0)/log(10.0d0), na)
	        agrid = agrid + ALOW - 1.0d0
	        kgrid = logspace(log(ALOW - 1.0d0*ALOW + 1.0d0)/log(10.0d0), log(khigh - 1.0d0*ALOW + 1.0d0)/log(10.0d0), nk)
	        kgrid = kgrid + ALOW - 1.0d0
        end if

		! setup spline
		invTa = spbas(ra,agrid)
		call dgetrf(ra,ra,invTa,ra,IPIV,INFO)
	    call dgetri(ra,invTa,ra,IPIV,WORK,ra,INFO)


        ! ********************************************************************************
		! INNER LOOP: value function iteration
        ! ********************************************************************************
		print *, " "
        distance = 1d+4
        howard = 0
        howardcnt = 0

        apmat = 0.0d0
	    awpmat = 0.0d0
	    anpmat = 0.0d0
	    lmat = 0.0d0

        call system_clock(ct1, cr)

        ! iteration = 1
        ! do while (distance > tol)
		do iteration = 1,maxiter

            !$omp parallel do default(shared) private(ie,ib,ix,betanow,xnow,EV,je,jb,jx,cfa,ia,anow,awp,anp,TW,TN)
			do ie = 1,ne

				ib = idmat(ie,1)
				ix = idmat(ie,2)
				betanow = bgrid(ib)
                xnow = xgrid(ix)

                ! compute conditional expected value functions on grid points
                EV = 0.0d0
                do je = 1,ne

					jb = idmat(je,1)
					jx = idmat(je,2)
                    EV = EV + pbij(ib,jb)*pxij(ix,jx)*reshape(Vmat(:,je),(/na/))

                end do

				! fit spline
				cfa = spfit(invTa,EV,ra,agrid)

                do ia = 1,na

                    anow = agrid(ia)
                    awp = awpmat(ia,ie)
                    anp = anpmat(ia,ie)

					call calcwnval(xnow,anow,r,w,Yagg,atilde,betanow,ALOW,agrid,cfa,tol,howard,awp,anp,TW,TN)

					TWmat(ia,ie) = TW
					TNmat(ia,ie) = TN
                    awpmatnew(ia,ie) = awp
                    anpmatnew(ia,ie) = anp

                end do

            end do
            !$omp end parallel do

            ! computing TV
            do ia = 1,na

                do ie = 1,ne

                    if (TWmat(ia,ie) > TNmat(ia,ie)) then
                        TVmat(ia,ie) = TWmat(ia,ie)
                        apmatnew(ia,ie) = awpmatnew(ia,ie)
                        lmatnew(ia,ie) = 1
                    else
                        TVmat(ia,ie) = TNmat(ia,ie)
                        apmatnew(ia,ie) = anpmatnew(ia,ie)
                        lmatnew(ia,ie) = 0
                    end if

                end do

            end do


            ! wnmax = maxval(maxval(TWmat-TNmat,1),1)
            ! wnmin = minval(minval(TWmat-TNmat,1),1)

			!if (logdiff==1) then
            ! 	distance = maxval(maxval(abs((log(TVmat)-log(Vmat))),1),1)
			!else
            	distance = maxval(maxval(abs(TVmat-Vmat),1),1)
            !end if
            distancevec(iteration) = distance

            distpolmax(1) = maxval(maxval(abs(awpmatnew-awpmat),1),1)
		    distpolmean(1) = sum(sum(abs(awpmatnew-awpmat),1),1)
		    distpolmax(2) = maxval(maxval(abs(anpmatnew-anpmat),1),1)
		    distpolmean(2) = sum(sum(abs(anpmatnew-anpmat),1),1)
		    distpolmax(3) = maxval(maxval(abs(lmatnew-lmat),1),1)
		    distpolmean(3) = sum(sum(abs(lmatnew-lmat),1),1)

            Vmat = TVmat
            Wmat = TWmat
            Nmat = TNmat
            apmat = apmatnew
		    awpmat = awpmatnew
		    anpmat = anpmatnew
		    lmat = lmatnew

            ! diagnosis
            if (mod(iteration,diagnum)==0) then

				! write(*,"('  iteration ', I4, '   ||TV-V|| = ', F8.5, ' max(W-N) = ', F8.5, ' min(W-N) = ', F8.5)") &
				! iteration, distance, wnmax, wnmin
				write(*,"('  iteration ', I4, '   ||TV-V|| = ', F10.8)") &
				iteration, distance
                ! write(*,"('          ||Tawp-awp||   ||Tanp-anp||   ||Tlmat-lmat|| = ( ', F8.5, ', ', F8.5, ', ', F8.5, ')')") distpolmax
			    ! !write(*,"('          sum|Tawp-awp|  sum|Tanp-anp|  sum|Tlmat-lmat|= ( ', F8.5, ', ', F8.5, ', ', F8.5, ')')") distpolmean
			    ! write(*,"('          Howard''s improvement algorithm: ', I2)") howardcnt
			    ! print *, howard

            end if

			if ((distance<=tol) .or. (iteration==maxiter)) then

				print *, " "
				if (iteration<maxiter) then
					print *, "  Value function iteration converged:"
				else
					print *, "  Maximum number of iterations reached:"
				end if
				do i=max(iteration-2,1),iteration
					write(*,"('  iteration ', I4, '   ||TV-V|| = ', F10.8)") &
	            	i, distancevec(i)
				end do

				exit

            end if

            !if (howflag==1) then
            !
            !    if ((maxval(distpolmax,1)<=1d-3) .and. (howard==0)) then
            !        ! print *, "  Howard''s improvement algorithm starts"
            !        ! pause
            !        howard = 1
            !
            !    end if
            !
            !    if (howard==1) then
            !
            !        howardcnt = howardcnt + 1
            !        if (howardcnt>10) then
            !            howard = 0
            !            howardcnt = 0
            !        end if
            !
            !    end if
            !
            !end if

            if (howflag==1) then

			    if (howard==0) then
				    ! print *, "  Howard''s improvement algorithm starts"
				    ! pause
				    howard = 1
			    end if

			    if (howard==1) then

				    howardcnt = howardcnt + 1
				    if (howardcnt>iteration/10) then ! NOTE: iteration dependent howard updating
					    howard = 0
					    howardcnt = 0
				    end if

			    end if

		    end if

        end do

        call system_clock(ct2, cr)
        eptime = dble(ct2-ct1)/cr
        write(*,"('  VFI Elasped time = ', F10.5)") eptime


        ! ********************************************************************************
        ! threshold asset for work decision
		print *, " "
		print *, " calculating threshold asset"
        minasset = 0.0d0

        do ie = 1,ne

            ib = idmat(ie,1)
            ix = idmat(ie,2)
			betanow = bgrid(ib)
            xnow = xgrid(ix)

            ! compute conditional expected value functions on grid points
            EV = 0.0d0
            do je = 1,ne
                jb = idmat(je,1)
                jx = idmat(je,2)
                EV = EV + pbij(ib,jb)*pxij(ix,jx)*reshape(Vmat(:,je),(/na/))
            end do

            ! fit spline
			cfa = spfit(invTa,EV,ra,agrid)

			minasset(ie) = bisectma(agrid,cfa,Wmat(:,ie),Nmat(:,ie),xnow,r,w,Yagg,atilde,betanow,ALOW,tol,1)

        end do ! ie loop

		! print *, " done."


        ! ********************************************************************************
        ! stationary distribution based on minasset
		if (transmat==0) then
		! iterarive method without transition matrix
        ! it is faster as we will not visit places with no population???

	        distance = 1d+4
	        ! iteration = 1
	        call system_clock(ct1, cr)

	        ! do while (distance > toldist)
			do iteration = 1,2*maxiter

	            !$omp parallel do private(ie,ib,ix,xnow,betanow,EV,je,jb,jx,cfa,ik,anow,Trw,Trn,kwhigh,knhigh,jk,wk,ap)
				do ie = 1,ne

	                ib = idmat(ie,1)
	                ix = idmat(ie,2)
				    betanow = bgrid(ib)
	                xnow = xgrid(ix)

	                ! compute conditional expected value functions on grid points
	                EV = 0.0d0
	                do je = 1,ne
	                    jb = idmat(je,1)
	                    jx = idmat(je,2)
	                    EV = EV + pbij(ib,jb)*pxij(ix,jx)*reshape(Vmat(:,je),(/na/))
	                end do

	                ! fit spline
				    cfa = spfit(invTa,EV,ra,agrid)

	                do ik = 1,nk

	                    anow = kgrid(ik)
	                    if (anow >= minasset(ie)) then
	                        ! Not work
							ndist(ik,ie) = 0.0d0
							call transfers(xnow,anow,r,w,0.0d0,Yagg,Trscale,Transprog,T0,Trn)
	                        knhigh = (1.0d0+r)*anow + Trn
	                        call gss(ALOW,knhigh,ap,tol**2,0,anow,xnow,agrid,cfa,w,r,Trn,betanow)
	                    else
	                        ! Work
							ndist(ik,ie) = HBAR
							call transfers(xnow,anow,r,w,HBAR,Yagg,Trscale,Transprog,T0,Trw)
	                        kwhigh = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Trw
	                        call gss(ALOW,kwhigh,ap,tol**2,1,anow,xnow,agrid,cfa,w,r,Trw,betanow)
	                    end if

	                    ! NOTE: finding indices and weights on grid points of k in the next periods
						! Ref. Young (2010?) and Rios-Rull (1999?)
						call findjw(ap,kgrid,jk,wk)
						jkmat(ik,ie) = jk
						wkmat(ik,ie) = wk

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

				if (logdiff==1) then
	            	distance = maxval(maxval(log(Tmumat)-log(mumat),1),1)
				else
	            	distance = maxval(maxval(Tmumat-mumat,1),1)
				end if
	            mumat = Tmumat
	            musum = sum(sum(mumat,1),1)

	            if (mod(iteration,diagnum)==0) then

					write(*,"('  iteration ', I4, '   ||Tmu-mu|| = ', F10.8, '  sum(mu) = ', F10.8)") &
					iteration, distance, musum

	            end if

				distancevec(iteration) = distance

				if ((distance<=toldist) .or. (iteration==maxiter)) then

					print *, " "
					if (iteration<maxiter) then
						print *, "  Distribution iteration converged:"
					else
						print *, "  Maximum number of iterations reached:"
					end if
					do i=max(iteration-2,1),iteration
						write(*,"('  iteration ', I4, '   ||Tmu-mu|| = ', F10.8)") &
						i, distancevec(i)
					end do

					exit

				end if
	            ! iteration = iteration + 1
	            ! if (iteration>2000) exit

	        end do

		else

			! with transition matrix
			!$omp parallel do private(ie,ib,ix,xnow,betanow,EV,je,jb,jx,cfa,ik,anow,Trw,Trn,kwhigh,knhigh,jk,wk,ap)
			do ie = 1,ne

				ib = idmat(ie,1)
				ix = idmat(ie,2)
				betanow = bgrid(ib)
				xnow = xgrid(ix)

				! compute conditional expected value functions on grid points
				EV = 0.0d0
				do je = 1,ne
					jb = idmat(je,1)
					jx = idmat(je,2)
					EV = EV + pbij(ib,jb)*pxij(ix,jx)*reshape(Vmat(:,je),(/na/))
				end do

				! fit spline
				cfa = spfit(invTa,EV,ra,agrid)

				do ik = 1,nk

                    anow = kgrid(ik)
                    if (anow >= minasset(ie)) then
                        ! Not work
						ndist(ik,ie) = 0.0d0
						call transfers(xnow,anow,r,w,0.0d0,Yagg,Trscale,Transprog,T0,Trn)
                        knhigh = (1.0d0+r)*anow + Trn
                        call gss(ALOW,knhigh,ap,tol**2,0,anow,xnow,agrid,cfa,w,r,Trn,betanow)
                    else
                        ! Work
						ndist(ik,ie) = HBAR
						call transfers(xnow,anow,r,w,HBAR,Yagg,Trscale,Transprog,T0,Trw)
                        kwhigh = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Trw
                        call gss(ALOW,kwhigh,ap,tol**2,1,anow,xnow,agrid,cfa,w,r,Trw,betanow)
                    end if

                    ! NOTE: finding indices and weights on grid points of k in the next periods
					! Ref. Young (2010?) and Rios-Rull (1999?)
					call findjw(ap,kgrid,jk,wk)
					jkmat(ik,ie) = jk
					wkmat(ik,ie) = wk

                end do

            end do
			!$omp end parallel do

			AA = 0.0d0
	        if (spblas==1) index = 0

	        do ie = 1,ne

	            ib = idmat(ie,1)
	            ix = idmat(ie,2)

	            do ik = 1,nk

	                if (spblas==1) csrRowPtr(nk*(ie-1)+ik) = index+1

	                do je = 1,ne

						jb = idmat(je,1)
			            jx = idmat(je,2)

                        AA(nk*(ie-1)+ik,nk*(je-1)+jkmat(ik,ie)) = pbij(ib,jb)*pxij(ix,jx)*wkmat(ik,ie)
						AA(nk*(ie-1)+ik,nk*(je-1)+jkmat(ik,ie)+1) = pbij(ib,jb)*pxij(ix,jx)*(1.0d0-wkmat(ik,ie))

						if (spblas==1) then

	                        index = index+1
	                        csrVal(index) = pbij(ib,jb)*pxij(ix,jx)*wkmat(ik,ie)
	                        csrColInd(index) = nk*(je-1)+jkmat(ik,ie)

		                    index = index+1
		                    csrVal(index) = pbij(ib,jb)*pxij(ix,jx)*(1.0d0-wkmat(ik,ie))
		                    csrColInd(index) = nk*(je-1)+jkmat(ik,ie)+1

						end if

	                end do

	            end do

	        end do

			if (spblas==1) then

		        csrRowPtr(nk*ne+1) = index+1
		        !   Create CSR matrix
		        i = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,nk*ne,nk*ne,csrRowPtr,csrRowPtr(2),csrColInd,csrVal)
		        !   Create matrix descriptor
		        descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
		        !   Analyze sparse matrix; chose proper kernels and workload balancing strategy
		        info = MKL_SPARSE_OPTIMIZE(csrA)

			end if

			! Option 1: iterative method
	        if (transmat==1) then

	            distance = 1d+4
	            iteration = 0

	            do while(distance>toldist) ! .and. iter<maxiter)

					if (spblas==1) then
	                	info = MKL_SPARSE_D_MV(SPARSE_OPERATION_TRANSPOSE,1.0d0,csrA,descrA,reshape(mumat,(/nk*ne/)),0.0d0,muvec)
					else
						! muvec = matmul(transpose(AA),reshape(mumat,(/nk*ne/)))
						call dgemv('T', nk*ne, nk*ne, 1.0d0, AA, nk*ne, reshape(mumat,(/nk*ne/)), 1, 0.0d0, muvec, 1)
					end if

	                Tmumat = reshape(muvec,(/nk,ne/))

					if (logdiff==1) then
		            	distance = maxval(maxval(log(Tmumat)-log(mumat),1),1)
					else
		            	distance = maxval(maxval(Tmumat-mumat,1),1)
					end if

		            mumat = Tmumat
					iteration = iteration + 1
		            musum = sum(sum(mumat,1),1)

		            if (mod(iteration,diagnum)==0) then

						write(*,"('  iteration ', I4, '   ||Tmu-mu|| = ', F10.8, '  sum(mu) = ', F10.8)") &
						iteration, distance, musum

		            end if

	            end do

	        ! Option 2: lapack (dgeev)
	        else if (transmat==2) then

	            call eig(transpose(AA),mumat)

	        ! ! Option 3: arpack (dnaupd and dneupd)
			! NOTE: need to install arpack library, see e.g. https://ylqk9.wordpress.com/build-arpack-for-linux-and-windows/
	        else if (transmat==3) then

	            call eigs(AA,csrA,descrA,1,'LM',toldist,mumat)

	        end if

	        !   Release internal representation of CSR matrix
	        info = MKL_SPARSE_DESTROY(csrA)

		end if ! transmat

        call system_clock(ct2, cr)
        eptime = dble(ct2-ct1)/cr
        write(*,"('  Distribution Elasped time = ', F10.5)") eptime


		! ****************************************
		! calculate aggregate variables
        muk = sum(mumat,2)

        ! NOTE: Ks, Ls and Ys are capital, labor and output implied by the supply side
        Ks = 0.0d0
        Ls = 0.0d0
        do ie = 1,ne
			! pick up the column in which idmat(ie,2) = ix
			Ks = Ks + sum(mumat(:,ie)*kgrid)
            Ls = Ls + sum(mumat(:,ie)*ndist(:,ie)*xgrid(idmat(ie,2)),1)
        end do

        Ys = zagg * (Ks**(ALPHA)) * (Ls**(1.0d0-ALPHA))
        EPratio = sum(sum(mumat*ndist/HBAR,1),1)

        Tragg = 0.0d0

        do ie = 1,ne

            ix = idmat(ie,2)
            xnow = xgrid(ix)

            do ik = 1,nk

                anow = kgrid(ik)

				! NOTE: Call tax and transfer functions. Progressive tax is to be done?
                !if (ndist(ik,ie)>0.0d0) then
                !    earns = w*xnow*HBAR
                !else
                !    earns = 0.0d0
                !end if

                !call taxfunc(earns,anow,r,Yagg,taul,a0,a1,a2,tax)
                call transfers(xnow,anow,r,w,ndist(ik,ie),Yagg,Trscale,Transprog,T0,Tr)
                !Taxagg = Taxagg + tax * mumat(ik,ie)
                Tragg = Tragg + Tr*mumat(ik,ie)

            end do

        end do

        ! Taxagg = taul*w*Ls
        ! KY = Ks/Ys
        ! TrY = Tragg/Ys ! Transfer/GDP ratio
        ! Govt= Taxagg-Tragg ! Total tax revenue - Total Transfers, Government expenditure is like a residual; should be positive

		!distanceout = max( abs(log(Ls)-log()), abs(log(Ks)-log(Kagg)) )
        !distanceout = max( abs(log(Ls)-log()), abs(log(Ks)-log(Kagg)) )
        ! we minimize K/L deviation
        distancedamp = log(Ks/Ls) - log(Kagg/Lagg)
        distanceout = abs(distancedamp)

		write(*,*) " "
		write(*,"('  ======> pct(Ks/Ls-Kagg/Lagg) = ', F10.8)") distanceout
		write(*,*) " "
		! write(*,"('  Aggregated statistics')")
		! write(*,"('   L   = ', F8.5, '  Ls = ', F8.5, ' K    = ', F8.5, ' Ks    =', F8.5)") Lagg, Ls, Kagg, Ks
		! write(*,"('   EP  = ', F8.5, ' K/Y = ', F8.5, ' Govt = ', F8.5)") EPratio, KY, Govt
        ! write(*,"('   Tax = ', F8.5, '  Tr = ', F8.5, ' Tr/Y = ', F8.5, ' Tax/Y =', F8.5)") Taxagg, Tragg, Tragg/Ys, Taxagg/Ys

		if (calcdist) then
			call calcdiststat(xgrid,kgrid,ndist,mumat,muk,idmat,w,r,HBAR,Trscale,Transprog,T0,EPratio,Yagg,Tragg, &
				Giniwage,Giniearns,Giniwealth,ShareW51,ShareW52,ShareW53,ShareW54,ShareW55, &
				EmpW51,EmpW52,EmpW53,EmpW54,EmpW55,TrInc51,TrInc52,TrInc53,TrInc54,TrInc55,TrW51,TrW52,TrW53,TrW54,TrW55)
			write(*,*) " "
			write(*,"('  Disaggregated statistics')")
			write(*,"('   Giniwage = ', F8.5, ' Giniearns = ', F8.5, ' Giniwealth = ', F8.5)") Giniwage, Giniearns, Giniwealth
	        write(*,"('   ShareW5-1st = ', F8.3, ' 2nd = ', F8.3, ' 3rd = ', F8.3, ' 4th = ', F8.3, ' 5th = ', F8.3)") ShareW51, ShareW52, ShareW53, ShareW54, ShareW55
	        write(*,"('   ErateW5-1st = ', F8.5, ' 2nd = ', F8.5, ' 3rd = ', F8.5, ' 4th = ', F8.5, ' 5th = ', F8.5)") EmpW51, EmpW52, EmpW53, EmpW54, EmpW55
            write(*,"('   TrW5-1st    = ', F8.3, ' 2nd = ', F8.3, ' 3rd = ', F8.3, ' 4th = ', F8.3, ' 5th = ', F8.3)") TrW51, TrW52, TrW53, TrW54, TrW55
	        write(*,"('   TrInc5-1st  = ', F8.3, ' 2nd = ', F8.3, ' 3rd = ', F8.3, ' 4th = ', F8.3, ' 5th = ', F8.3)") TrInc51, TrInc52, TrInc53, TrInc54, TrInc55
		end if

		!write(*,*) " "
		!write(*,"('  dampling weight =', F10.5)") damp
		!print *,'===================================================================='

        ! updating price
        ! if Ks/Ls is too high, r should go down
        ! if Ks/Ls is too low, r should increase
        ! r = (1.0d0 - (1.0d0-damp)*distancedamp)*r

        Kagg = Ks
        Lagg = Ls
        rnew = ALPHA*zagg*(Kagg/Lagg)**(ALPHA-1.0d0) - DELTA
        wnew = (1.0d0-ALPHA)*zagg*(Kagg/Lagg)**ALPHA
        r = damp*r + (1.0d0-damp)*rnew
        w = damp*w + (1.0d0-damp)*wnew

        iterationout = iterationout + 1
		! if (iterationout > 100) exit

    end do


    do ie = 1,ne
        write(*,"('  minasset for xgrid ', I4, '   is = ', F8.2)") ie, minasset(ie)
    end do
	call system_clock(cto2, cro)
	eptimeo = dble(cto2-cto1)/cro
	write(*,"('  Total Elasped time = ', F10.5)") eptimeo


end subroutine calcss


function bisectma(agrid,cfa,Wvec,Nvec,xnow,r,w,Yagg,atilde,betanow,ALOW,tol,ascend) result(ma0)


	real(8), intent(in) :: agrid(:), cfa(:,:), Wvec(:), Nvec(:), xnow, r, w, Yagg, atilde, betanow, ALOW, tol
	integer, intent(in) :: ascend
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

			! do while (distance > tol**2.0d0) ! NOTE 180507 old calibration does not converge when tol*10.0d0
	        do while (distance > tol*10.0d0)

	            ma0 = (malow+mahigh)/2.0d0
	            ! anow = ma0

				! call calcwnval(xnow,ma0,r,w,Yagg,atilde,HBAR,Trscale,Transprog,T0,betanow,B0,taul,ALOW,agrid,cfa,tol,0,0,awp,anp,TW,TN)
				call calcwnval(xnow,ma0,r,w,Yagg,atilde,betanow,ALOW,agrid,cfa,tol,0,awp,anp,TW,TN)

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


subroutine calcwnval(xnow,anow,r,w,Yagg,atilde,betanow,ALOW,agrid,cfa,tol,howard,awp,anp,TW,TN)
      ! call calcwnval(xnow,anow,r,w,Yagg,atilde,betanow,ALOW,agrid,cfa,tol,howard,awp,anp,TW,TN)
! call calcwnval(xnow,anow,r,w,Yagg,atilde,betanow,agrid,cfa,tol,0,howard,awp,anp,TW,TN)

	real(8), intent(in) :: xnow, anow, r, w, Yagg, atilde, betanow, ALOW
	real(8), intent(in) :: agrid(:), cfa(:,:), tol
	integer, intent(in) :: howard
	real(8), intent(inout) :: awp, anp
	real(8), intent(out) :: TW, TN
	real(8) Trw, kwhigh, Trn, knhigh


	! value of work
	call transfers(xnow,anow,r,w,HBAR,Yagg,Trscale,Transprog,T0,Trw)
	kwhigh = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Trw
	if ((1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow - alow <= 0.0d0) then
		write(*,"('  max consumption (when work) is negative', F10.2)") kwhigh
	end if

	! NOTE: choose cubic spline or linear interpolation to be done
	! option1: cubic spline
	if (howard==0) call gss(ALOW,kwhigh,awp,tol**2,1,anow,xnow,agrid,cfa,w,r,Trw,betanow)
	TW = -objneg_w(awp,anow,xnow,agrid,cfa,w,r,Trw,betanow)
	! option2: linear interpolation
	!call gss(ALOW,kwhigh,awpmat(ia,ie),tol**2,1,anow,xnow,agrid,EV,w,r,Tr,betanow,B0,Hbar,taul,Trw)
	!TWmat(ia,ie) = -objneg_w(awpmat(ia,ie),anow,xnow,agrid,EV,w,r,Tr,betanow,B0,Hbar,taul,Trw)

	! value of non-employed
	if (anow>atilde) then
		call transfers(xnow,anow,r,w,0.0d0,Yagg,Trscale,Transprog,T0,Trn)
		knhigh = (1.0d0+r)*anow + Trn
		! option1: cubic spline
		if (howard==0) call gss(ALOW,knhigh,anp,tol**2,0,anow,xnow,agrid,cfa,w,r,Trn,betanow)
		TN = -objneg_n(anp,anow,agrid,cfa,r,Trn,betanow)
		! option2: linear interpolation
		!call gss(ALOW,knhigh,anpmat(ia,ie),tol**2,0,anow,xnow,agrid,EV,w,r,Trn,betanow,B0,Hbar,taul)
		!TNmat(ia,ie) = -objneg_n(anpmat(ia,ie),anow,agrid,EV,r,Trn,betanow)
	else
		anp = alow
		TN = -1d+20
	end if


end subroutine calcwnval


subroutine gss(xmin,xmax,x,critg,flag,anow,xnow,agrid,cfa,w,r,Tr,BETA)


	real(8), intent(in) :: xmin, xmax, critg, anow, xnow, w, r, Tr, beta
    real(8), intent(in) :: agrid(:)
    real(8), intent(in) :: cfa(:,:)
    ! for linear interpolation
    !real(8), intent(in) :: cfa(:)
    integer, intent(in) :: flag
	real(8), intent(out) :: x

	real(8) rg, a, b, c, d, z, fc, fd, crit, diff, fa, fb
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    if (flag==1) then
        fc = objneg_w(c,anow,xnow,agrid,cfa,w,r,Tr,BETA)
        fa = objneg_w(a,anow,xnow,agrid,cfa,w,r,Tr,BETA)
        fb = objneg_w(b,anow,xnow,agrid,cfa,w,r,Tr,BETA)
    else
        fc = objneg_n(c,anow,agrid,cfa,r,Tr,BETA)
        fa = objneg_n(a,anow,agrid,cfa,r,Tr,BETA)
        fb = objneg_n(b,anow,agrid,cfa,r,Tr,BETA)
    end if

    ! x = c

	d = a + (1.0d0-rg)*(b-a)
    if (flag==1) then
        fd = objneg_w(d,anow,xnow,agrid,cfa,w,r,Tr,beta)
    else
        fd = objneg_n(d,anow,agrid,cfa,r,Tr,beta)
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
            if (flag==1) then
                fd = objneg_w(d,anow,xnow,agrid,cfa,w,r,Tr,beta)
            else
                fd = objneg_n(d,anow,agrid,cfa,r,Tr,beta)
            end if

	    else

        	z = a + rg*(d-a)
        	b = d
        	d = c
        	fd = fc
        	c = z
            if (flag==1) then
                fc = objneg_w(c,anow,xnow,agrid,cfa,w,r,Tr,beta)
            else
                fc = objneg_n(c,anow,agrid,cfa,r,Tr,beta)
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


function objneg_w(ap,anow,xnow,agrid,cfa,w,r,Tr,BETA) result(f)


    real(8), intent(in) :: ap, anow, xnow, w, r, Tr, beta !, B0, HBAR, taul
    real(8), intent(in) :: agrid(:)
    real(8), intent(in) :: cfa(:,:)
    ! for linear interpolation
    !real(8), intent(in) :: cfa(:)
    real(8) cnow, ev, edv, ed2v, f
	integer ra

	ra = size(agrid,1)-2


    cnow = (1.0d0-taul)*w*xnow*HBAR + (1.0d0+r)*anow + Tr - ap

    ! option1: cubic spline
    call speva(cfa,ap,ra,agrid,EV,EDV,ED2V)
	f = log(cnow) - B0 + beta*EV

    ! option2: linear interpolation
    !f = log(cnow) - B0 + beta*lini(agrid,cfa,ap)

    f = -f


end function objneg_w


function objneg_n(ap,anow,agrid,cfa,r,Tr,beta) result(f)


    real(8), intent(in) :: ap, anow, r, Tr, beta
    real(8), intent(in) :: agrid(:)
    real(8), intent(in) :: cfa(:,:)
    ! for linear interpolation
    !real(8), intent(in) :: cfa(:)
    real(8) cnow, ev, edv, ed2v, f
	integer ra

	ra = size(agrid,1)-2

    cnow = (1.0d0+r)*anow + Tr - ap

    ! option1: cubic spline
    call speva(cfa,ap,ra,agrid,EV,EDV,ED2V)
	f = log(cnow) + beta*EV

    ! option2: linear interpolation
	!f = log(cnow) + beta*lini(agrid,cfa,ap)

    f = -f


end function objneg_n


subroutine calcdiststat(xgrid,kgrid,ndist,mumat,muk,idmat,w,r,HBAR,Trscale,Transprog,T0,EPratio,Yagg,Tragg, &
	Giniwage,Giniearns,Giniwealth,ShareW51,ShareW52,ShareW53,ShareW54,ShareW55, &
	EmpW51,EmpW52,EmpW53,EmpW54,EmpW55,TrInc51,TrInc52,TrInc53,TrInc54,TrInc55,TrW51,TrW52,TrW53,TrW54,TrW55)


	real(8), intent(in) :: xgrid(:), kgrid(:), ndist(:,:), mumat(:,:), muk(:)
	integer, intent(in) :: idmat(:,:)
	real(8), intent(in) :: w, r, HBAR, EPratio, Yagg, Tragg, Trscale, Transprog, T0
	real(8), intent(out) :: Giniwage, Giniearns, Giniwealth, ShareW51, ShareW52, ShareW53, ShareW54, ShareW55, &
		EmpW51, EmpW52, EmpW53, EmpW54, EmpW55, TrInc51, TrInc52, TrInc53, TrInc54, TrInc55,TrW51,TrW52,TrW53,TrW54,TrW55
	integer ne, nk, ie, ix, ik, ii
	real(8), allocatable :: incgrid(:), muinc(:), trgrid(:)
	real(8) W5threshold(4), cummmk, wgtik(4), negW, SumwgtW5(5), EsumW5(5), ShareW5(5), ShareE5(5), TrsumW5(5), trans !, ShareW10(10), ShareW101
    real(8) Inc5threshold(4), SumwgtInc5(5), TrsumInc5(5)


	ne = size(idmat,1)
	nk = size(kgrid,1)
	allocate(incgrid(nk*ne),muinc(nk*ne),trgrid(nk*ne))

	! distribution (maybe optional)
	call calwagegini(xgrid,ndist,mumat,idmat,w,HBAR,EPratio,Giniwage)
	call calearnsgini(xgrid,ndist,mumat,idmat,w,Giniearns)
	!call calwealthgini(kgrid,mumat,Giniwealth)
	call ginicoef(muk,kgrid,nk,5,Giniwealth,ShareW5,negW)
	! call ginicoef(muk,kgrid,nk,10,Giniwealth,ShareW10,negW)

	! ShareW101 = ShareW10(1)*100.0d0

	! Employment rate by wealth
	! calculating threshold
	W5threshold = 0         ! index
	cummmk = 0.0d0
	wgtik = 0.0d0           ! calculating division rule for threshold

	do ik = 1,nk
		cummmk = cummmk + muk(ik)
		if (cummmk >= 0.2d0 .and. W5threshold(1) == 0 ) then
			W5threshold(1) = ik
			wgtik(1) = 1.0d0 - (cummmk - 0.2d0)/(muk(ik))
		elseif (cummmk >= 0.4d0 .and. W5threshold(2) == 0 ) then
			W5threshold(2) = ik
			wgtik(2) = 1.0d0 - (cummmk - 0.4d0)/(muk(ik))
		elseif (cummmk >= 0.6d0 .and. W5threshold(3) == 0 ) then
			W5threshold(3) = ik
			wgtik(3) = 1.0d0 - (cummmk - 0.6d0)/(muk(ik))
		elseif (cummmk >= 0.8d0 .and. W5threshold(4) == 0 ) then
			W5threshold(4) = ik
			wgtik(4) = 1.0d0 - (cummmk - 0.8d0)/(muk(ik))
		end if
	end do

	SumwgtW5 = 0.0d0
	EsumW5 = 0.0d0
    TrsumW5 = 0.0d0
	ShareW5 = 0.0d0

	do ik = 1,nk
		do ie = 1,ne
            ix = idmat(ie,2)
            call transfers(xgrid(ix),kgrid(ik),r,w,ndist(ik,ie),Yagg,Trscale,Transprog,T0,trans)

			if (ik < W5threshold(1)) then
				SumwgtW5(1) = SumwgtW5(1) + mumat(ik,ie)
				EsumW5(1) = EsumW5(1) + (ndist(ik,ie)/hbar)*mumat(ik,ie)
				ShareW5(1) = ShareW5(1) + kgrid(ik)*mumat(ik,ie)
                TrsumW5(1) = TrsumW5(1) + trans*mumat(ik,ie)
			elseif (ik == W5threshold(1)) then
				SumwgtW5(1) = SumwgtW5(1) + mumat(ik,ie)*wgtik(1)
				SumwgtW5(2) = SumwgtW5(2) + mumat(ik,ie)*(1.0d0-wgtik(1))
				EsumW5(1) = EsumW5(1) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*wgtik(1)
				EsumW5(2) = EsumW5(2) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*(1.0d0-wgtik(1))
				ShareW5(1) = ShareW5(1) + kgrid(ik)*mumat(ik,ie)*wgtik(1)
				ShareW5(2) = ShareW5(2) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(1))
                TrsumW5(1) = TrsumW5(1) + trans*mumat(ik,ie)*wgtik(1)
                TrsumW5(2) = TrsumW5(2) + trans*mumat(ik,ie)*(1.0d0-wgtik(1))
			elseif (ik > W5threshold(1) .and. ik < W5threshold(2)) then
				SumwgtW5(2) = SumwgtW5(2) + mumat(ik,ie)
				EsumW5(2) = EsumW5(2) + (ndist(ik,ie)/hbar)*mumat(ik,ie)
				ShareW5(2) = ShareW5(2) + kgrid(ik)*mumat(ik,ie)
                TrsumW5(2) = TrsumW5(2) + trans*mumat(ik,ie)
			elseif (ik == W5threshold(2)) then
				SumwgtW5(2) = SumwgtW5(2) + mumat(ik,ie)*wgtik(2)
				SumwgtW5(3) = SumwgtW5(3) + mumat(ik,ie)*(1.0d0-wgtik(2))
				EsumW5(2) = EsumW5(2) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*wgtik(2)
				EsumW5(3) = EsumW5(3) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*(1.0d0-wgtik(2))
				ShareW5(2) = ShareW5(2) + kgrid(ik)*mumat(ik,ie)*wgtik(2)
				ShareW5(3) = ShareW5(3) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(2))
                TrsumW5(2) = TrsumW5(2) + trans*mumat(ik,ie)*wgtik(2)
                TrsumW5(3) = TrsumW5(3) + trans*mumat(ik,ie)*(1.0d0-wgtik(2))
			elseif (ik > W5threshold(2) .and. ik < W5threshold(3)) then
				SumwgtW5(3) = SumwgtW5(3) + mumat(ik,ie)
				EsumW5(3) = EsumW5(3) + (ndist(ik,ie)/hbar)*mumat(ik,ie)
				ShareW5(3) = ShareW5(3) + kgrid(ik)*mumat(ik,ie)
                TrsumW5(3) = TrsumW5(3) + trans*mumat(ik,ie)
			elseif (ik == W5threshold(3)) then
				SumwgtW5(3) = SumwgtW5(3) + mumat(ik,ie)*wgtik(3)
				SumwgtW5(4) = SumwgtW5(4) + mumat(ik,ie)*(1.0d0-wgtik(3))
				EsumW5(3) = EsumW5(3) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*wgtik(3)
				EsumW5(4) = EsumW5(4) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*(1.0d0-wgtik(3))
				ShareW5(3) = ShareW5(3) + kgrid(ik)*mumat(ik,ie)*wgtik(3)
				ShareW5(4) = ShareW5(4) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(3))
                TrsumW5(3) = TrsumW5(3) + trans*mumat(ik,ie)*wgtik(3)
                TrsumW5(4) = TrsumW5(4) + trans*mumat(ik,ie)*(1.0d0-wgtik(3))
			elseif (ik > W5threshold(3) .and. ik < W5threshold(4)) then
				SumwgtW5(4) = SumwgtW5(4) + mumat(ik,ie)
				EsumW5(4) = EsumW5(4) + (ndist(ik,ie)/hbar)*mumat(ik,ie)
				ShareW5(4) = ShareW5(4) + kgrid(ik)*mumat(ik,ie)
                TrsumW5(4) = TrsumW5(4) + trans*mumat(ik,ie)
			elseif (ik == W5threshold(4)) then
				SumwgtW5(4) = SumwgtW5(4) + mumat(ik,ie)*wgtik(4)
				SumwgtW5(5) = SumwgtW5(5) + mumat(ik,ie)*(1.0d0-wgtik(4))
				EsumW5(4) = EsumW5(4) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*wgtik(4)
				EsumW5(5) = EsumW5(5) + (ndist(ik,ie)/hbar)*mumat(ik,ie)*(1.0d0-wgtik(4))
				ShareW5(4) = ShareW5(4) + kgrid(ik)*mumat(ik,ie)*wgtik(4)
				ShareW5(5) = ShareW5(5) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(4))
                TrsumW5(4) = TrsumW5(4) + trans*mumat(ik,ie)*wgtik(4)
                TrsumW5(5) = TrsumW5(5) + trans*mumat(ik,ie)*(1.0d0-wgtik(4))
			else
				SumwgtW5(5) = SumwgtW5(5) + mumat(ik,ie)
				EsumW5(5) = EsumW5(5) + (ndist(ik,ie)/hbar)*mumat(ik,ie)
				ShareW5(5) = ShareW5(5) + kgrid(ik)*mumat(ik,ie)
                TrsumW5(5) = TrsumW5(5) + trans*mumat(ik,ie)
			end if
		end do
	end do

	EmpW51 = EsumW5(1)/SumwgtW5(1)
	EmpW52 = EsumW5(2)/SumwgtW5(2)
	EmpW53 = EsumW5(3)/SumwgtW5(3)
	EmpW54 = EsumW5(4)/SumwgtW5(4)
	EmpW55 = EsumW5(5)/SumwgtW5(5)

	ShareW51 = 100.0d0*ShareW5(1)/sum(ShareW5)
	ShareW52 = 100.0d0*ShareW5(2)/sum(ShareW5)
	ShareW53 = 100.0d0*ShareW5(3)/sum(ShareW5)
	ShareW54 = 100.0d0*ShareW5(4)/sum(ShareW5)
	ShareW55 = 100.0d0*ShareW5(5)/sum(ShareW5)

    TrW51 = (TrsumW5(1)/SumwgtW5(1))/(Tragg)
	TrW52 = (TrsumW5(2)/SumwgtW5(2))/(Tragg)
	TrW53 = (TrsumW5(3)/SumwgtW5(3))/(Tragg)
	TrW54 = (TrsumW5(4)/SumwgtW5(4))/(Tragg)
	TrW55 = (TrsumW5(5)/SumwgtW5(5))/(Tragg)

	! Transfers by income
	! construct income distribution (not sorted)
	muinc = 0.0d0
	incgrid = 0.0d0
	trgrid = 0.0d0

	do ik = 1,nk
		do ie = 1,ne
			ii = ik + (ie-1)*nk
			muinc(ii) = mumat(ik,ie)
			ix = idmat(ie,2)
			if (kgrid(ik) > 0) then
				incgrid(ii) = w*xgrid(ix)*ndist(ik,ie) + r*kgrid(ik)
			else
				incgrid(ii) = w*xgrid(ix)*ndist(ik,ie)
			end if

			call transfers(xgrid(ix),kgrid(ik),r,w,ndist(ik,ie),Yagg,Trscale,Transprog,T0,trgrid(ii))
		end do
	end do

	! ordering income distirbution
	call sortdist((nk*ne), incgrid, muinc, trgrid, incgrid, muinc, trgrid)

	! calculating threshold
	Inc5threshold = 0         ! index
	cummmk = 0.0d0
	wgtik = 0.0d0           ! calculating division rule for threshold

	do ii = 1,(nk*ne)
		cummmk = cummmk + muinc(ii)
		if (cummmk >= 0.2d0 .and. Inc5threshold(1) == 0 ) then
			Inc5threshold(1) = ii
			wgtik(1) = 1.0d0 - (cummmk - 0.2d0)/(muinc(ii))
		elseif (cummmk >= 0.4d0 .and. Inc5threshold(2) == 0 ) then
			Inc5threshold(2) = ii
			wgtik(2) = 1.0d0 - (cummmk - 0.4d0)/(muinc(ii))
		elseif (cummmk >= 0.6d0 .and. Inc5threshold(3) == 0 ) then
			Inc5threshold(3) = ii
			wgtik(3) = 1.0d0 - (cummmk - 0.6d0)/(muinc(ii))
		elseif (cummmk >= 0.8d0 .and. Inc5threshold(4) == 0 ) then
			Inc5threshold(4) = ii
			wgtik(4) = 1.0d0 - (cummmk - 0.8d0)/(muinc(ii))
		end if
	end do

	SumwgtInc5 = 0.0d0
	TrsumInc5 = 0.0d0

	do ii = 1,(nk*ne)
		if (ii < Inc5threshold(1)) then
			SumwgtInc5(1) = SumwgtInc5(1) + muinc(ii)
			TrsumInc5(1) = TrsumInc5(1) + muinc(ii)*trgrid(ii)
		elseif (ii == Inc5threshold(1)) then
			SumwgtInc5(1) = SumwgtInc5(1) + muinc(ii)*wgtik(1)
			SumwgtInc5(2) = SumwgtInc5(2) + muinc(ii)*(1.0d0-wgtik(1))
			TrsumInc5(1) = TrsumInc5(1) + muinc(ii)*trgrid(ii)*wgtik(1)
			TrsumInc5(2) = TrsumInc5(2) + muinc(ii)*trgrid(ii)*(1.0d0-wgtik(1))
		elseif (ii > Inc5threshold(1) .and. ii < Inc5threshold(2)) then
			SumwgtInc5(2) = SumwgtInc5(2) + muinc(ii)
			TrsumInc5(2) = TrsumInc5(2) + muinc(ii)*trgrid(ii)
		elseif (ii == Inc5threshold(2)) then
			SumwgtInc5(2) = SumwgtInc5(2) + muinc(ii)*wgtik(2)
			SumwgtInc5(3) = SumwgtInc5(3) + muinc(ii)*(1.0d0-wgtik(2))
			TrsumInc5(2) = TrsumInc5(2) + muinc(ii)*trgrid(ii)*wgtik(2)
			TrsumInc5(3) = TrsumInc5(3) + muinc(ii)*trgrid(ii)*(1.0d0-wgtik(2))
		elseif (ii > Inc5threshold(2) .and. ii < Inc5threshold(3)) then
			SumwgtInc5(3) = SumwgtInc5(3) + muinc(ii)
			TrsumInc5(3) = TrsumInc5(3) + muinc(ii)*trgrid(ii)
		elseif (ii == Inc5threshold(3)) then
			SumwgtInc5(3) = SumwgtInc5(3) + muinc(ii)*wgtik(3)
			SumwgtInc5(4) = SumwgtInc5(4) + muinc(ii)*(1.0d0-wgtik(3))
			TrsumInc5(3) = TrsumInc5(3) + muinc(ii)*trgrid(ii)*wgtik(3)
			TrsumInc5(4) = TrsumInc5(4) + muinc(ii)*trgrid(ii)*(1.0d0-wgtik(3))
		elseif (ii > Inc5threshold(3) .and. ii < Inc5threshold(4)) then
			SumwgtInc5(4) = SumwgtInc5(4) + muinc(ii)
			TrsumInc5(4) = TrsumInc5(4) + muinc(ii)*trgrid(ii)
		elseif (ii == Inc5threshold(4)) then
			SumwgtInc5(4) = SumwgtInc5(4) + muinc(ii)*wgtik(4)
			SumwgtInc5(5) = SumwgtInc5(5) + muinc(ii)*(1.0d0-wgtik(4))
			TrsumInc5(4) = TrsumInc5(4) + muinc(ii)*trgrid(ii)*wgtik(4)
			TrsumInc5(5) = TrsumInc5(5) + muinc(ii)*trgrid(ii)*(1.0d0-wgtik(4))
		else
			SumwgtInc5(5) = SumwgtInc5(5) + muinc(ii)
			TrsumInc5(5) = TrsumInc5(5) + muinc(ii)*trgrid(ii)
		end if
	end do

	TrInc51 = (TrsumInc5(1)/SumwgtInc5(1))/(Tragg)
	TrInc52 = (TrsumInc5(2)/SumwgtInc5(2))/(Tragg)
	TrInc53 = (TrsumInc5(3)/SumwgtInc5(3))/(Tragg)
	TrInc54 = (TrsumInc5(4)/SumwgtInc5(4))/(Tragg)
	TrInc55 = (TrsumInc5(5)/SumwgtInc5(5))/(Tragg)



end subroutine calcdiststat


subroutine calwagegini(xgrid,ndist,mumat,idmat,wage,HBAR,EPratio,wagegini)

    implicit none

    real(8), intent(in) :: xgrid(:), ndist(:,:), mumat(:,:), wage, HBAR, EPratio
    integer, intent(in) :: idmat(:,:)
    real(8), intent(out) :: wagegini

    real(8), allocatable :: stnrpdf(:), stnrcdf(:), tempvec1(:), tempvec2(:)
    real(8), allocatable :: dc(:), db(:), tempmat1(:,:), tempmat2(:,:)
    integer :: ie, ix, nk, ne, nx

    nk = size(mumat,1)
    ne = size(mumat,2)
    nx = size(xgrid,1)
    allocate(stnrpdf(nx), stnrcdf(nx), tempvec1(nx), tempvec2(nx))
    allocate(dc(nx-1), db(nx-1), tempmat1(nk,ne), tempmat2(nk,ne))

    ! Calculate binary work choice (1 or 0) on the pdf
    tempmat1 = ndist/HBAR;
    ! Compute the conditional pdf of workers
    tempmat2 = (tempmat1*mumat)/EPratio

    ! stnrpdf: integrating capital out
    ! call calcstnrpdf(tempmat2,idmat,stnrpdf,stnrcdf)
    stnrpdf = 0.0d0
    do ie = 1,ne
        stnrpdf(idmat(ie,2)) = stnrpdf(idmat(ie,2)) + sum(tempmat2(:,ie))
    end do

    ! Generate the Lorenz Curve (tempvec2)
    tempvec1 = wage*xgrid*stnrpdf

    do ix = 1,nx
    	tempvec2(ix) = sum(tempvec1(1:ix))/sum(tempvec1)
    end do

    call calcgini(tempvec2,stnrpdf,wagegini)


end subroutine calwagegini


subroutine calearnsgini(xgrid,ndist,mumat,idmat,wage,earnsgini)

    implicit none

    real(8), intent(in) :: xgrid(:), ndist(:,:), mumat(:,:), wage
    integer, intent(in) :: idmat(:,:)
    real(8), intent(out) :: earnsgini

    real(8), allocatable :: stnrpdf(:), stnrcdf(:), tempvec1(:), tempvec2(:)
    real(8), allocatable :: dc(:), db(:), tempmat1(:,:), tempmat2(:,:)
    integer :: ie, ix, nk, ne, nx

    nk = size(mumat,1)
    ne = size(mumat,2)
    nx = size(xgrid,1)
    allocate(stnrpdf(nx), stnrcdf(nx), tempvec1(nx), tempvec2(nx))
    allocate(dc(nx-1), db(nx-1), tempmat1(nk,ne), tempmat2(nk,ne))

    tempmat2 = mumat

    ! stnrpdf: integrating capital out
    ! call calcstnrpdf(tempmat2,idmat,stnrpdf,stnrcdf)
    stnrpdf = 0.0d0
    do ie = 1,ne
        stnrpdf(idmat(ie,2)) = stnrpdf(idmat(ie,2)) + sum(tempmat2(:,ie))
    end do

    ! Generate the Lorenz Curve (tempvec2)
    tempmat1 = wage*ndist*mumat

    do ix = 1,nx
        tempvec1(ix) = sum(tempmat1(:,ix))
    end do

    do ix = 1,nx
    	tempvec2(ix) = sum(tempvec1(1:ix))/sum(tempvec1)
    end do

    call calcgini(tempvec2,stnrpdf,earnsgini)


end subroutine calearnsgini


subroutine calwealthgini(kgrid,mumat,wealthgini)

    implicit none

    real(8), intent(in) :: kgrid(:), mumat(:,:)
    real(8), intent(out) :: wealthgini

    real(8), allocatable :: stnrpdf(:), stnrcdf(:), tempvec1(:), tempvec2(:)
    real(8), allocatable :: dc(:), db(:), tempmat1(:,:), tempmat2(:,:)
    integer :: ie, ix, ik, nk, ne, nx

    nk = size(mumat,1)
    ne = size(mumat,2)
    ! nx = size(xgrid,1)
    allocate(stnrpdf(nk), stnrcdf(nk), tempvec1(nk), tempvec2(nk))
    allocate(dc(nk-1), db(nk-1), tempmat1(ne,nk), tempmat2(ne,nk))

    tempmat2 = transpose(mumat) ! ne x nk

    ! stnrpdf: integrating capital out
    stnrpdf = 0.0d0
    do ik = 1,nk
        stnrpdf(ik)=sum(tempmat2(:,ik))
    end do
    ! call calcstnrpdf(tempmat2,stnrpdf,stnrcdf)

    ! Generate the Lorenz Curve (tempvec2)
    tempvec1 = kgrid*stnrpdf
    ! aaa = sum(tempvec1)
    do ik = 1,nk
        tempvec2(ik) = sum(tempvec1(1:ik))/sum(tempvec1)
    end do

    call calcgini(tempvec2,stnrpdf,wealthgini)


end subroutine calwealthgini


subroutine calcgini(tempvec2,stnrpdf,gini)

    implicit none

    real(8), intent(in) :: tempvec2(:), stnrpdf(:)
    real(8), intent(out) :: gini
    integer ix, nx
    real(8), allocatable :: dc(:), db(:)
    real(8) aaa, bbb

    nx = size(tempvec2,1)
    allocate(dc(nx),db(nx))

    ! Calculate the gini coefficient of wage
    do ix = 1,nx-1
        dc(ix) = stnrpdf(ix+1) !stnrcdf(ix+1)-stnrcdf(ix)
        db(ix) = (tempvec2(ix+1)+tempvec2(ix))*0.5d0
    end do

    bbb = sum(dc*db)
    aaa = 0.5d0-bbb
    gini = 2.0d0*aaa


end subroutine calcgini


subroutine eig(AA,mu0)


    real(8), intent(in) :: AA(:,:)
    real(8), intent(out) :: mu0(:,:)
    integer info, lda, ldvr, lwork, n, nb, nx
    ! integer lda, ldvr, lwork
    integer, parameter :: nblock = 2048*64
    ! .. Local Arrays ..
    real(8), allocatable :: vr(:,:), wi(:), work(:), wr(:)
    real(8) dummy(1,1)


    nb = size(mu0,1)
    nx = size(mu0,2)
    n = nb*nx
    lda = n
    ldvr = n
    allocate(vr(ldvr,n),wi(n),wr(n))

	! print *, n
	! pause

    lwork = -1
    call dgeev('N', 'V', n, AA, lda, wr, wi, dummy, 1, vr, ldvr, dummy, lwork, info)

    ! Make sure that there is enough workspace for block size nb.
    lwork = max((nblock+2)*n,nint(dummy(1,1)))
    ! print *, lwork, nint(dummy(1,1))
    allocate(work(lwork))

    ! Compute the eigenvalues and right eigenvectors of A
    call dgeev('N', 'V', n, AA, lda, wr, wi, dummy, 1, vr, ldvr, work, lwork, info)
    ! print *, maxval(wr), sum(vr(:,maxloc(wr)))

    mu0 = reshape(vr(:,maxloc(wr))/sum(vr(:,maxloc(wr))),(/nb,nx/))


end subroutine eig


subroutine eigs(AA,csrA,descrA,nev,which,tol,mu0)

    ! for arpack
    real(8), intent(in) :: AA(:,:), tol
    !   Matrix descriptor
    TYPE(MATRIX_DESCR), intent(in) :: descrA     ! Sparse matrix descriptor
    !   CSR matrix representation
    TYPE(SPARSE_MATRIX_T), intent(in) :: csrA    ! Structure with sparse matrix
    integer, intent(in) :: nev
    character, intent(in) :: which*2
    real(8), intent(out) :: mu0(:,:)
    integer maxn, maxnev, maxncv, ldv
    parameter (maxn=100000, maxnev=4, maxncv=200, ldv=maxn)
    integer iparam(11), ipntr(14)
    logical select(maxncv)
    real(8) ax(maxn), d(maxncv,3), resid(maxn), v(ldv,maxncv), workd(3*maxn), workev(3*maxncv), workl(3*maxncv*maxncv+6*maxncv)
    integer ido, n, ncv, lworkl, info, ierr, j, ishfts, maxitr, mode1, nconv, converged, nb, nx
    real(8) sigmar, sigmai


    nb = size(mu0,1)
    nx = size(mu0,2)
    n = nb*nx
    ! parameters for arpack
    ! nev   = 1
    ncv = 100 !nev+2

    lworkl = 3*ncv**2+6*ncv
    ! tol    = 1d-5 !1d-10
    ido    = 0
    info   = 0

    ishfts = 1
    maxitr = 300
    mode1 = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1

    converged = 0

    do while (converged==0)

        call dnaupd(ido, 'I', n, which, nev, tol, resid, &
        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)

        if (ido .eq. -1 .or. ido .eq. 1) then

            ! workd(ipntr(2):ipntr(2)+n-1) = matmul(transpose(AA),workd(ipntr(1):ipntr(1)+n-1))
            ! call dgemv('T', n, n, 1.0d0, AA, n, workd(ipntr(1):ipntr(1)+n-1), 1, 0.0d0, workd(ipntr(2):ipntr(2)+n-1), 1)
            info = MKL_SPARSE_D_MV(SPARSE_OPERATION_TRANSPOSE,1.0d0,csrA,descrA,workd(ipntr(1):ipntr(1)+n-1),0.0d0,workd(ipntr(2):ipntr(2)+n-1))

        else if (ido==99) then

            converged = 1

        end if

    end do

    call dneupd(.true., 'A', select, d, d(1,2), v, ldv, sigmar, sigmai, workev, 'I', n, which, nev, tol, &
    resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr)

    mu0 = reshape(v(:,1)/sum(v(:,1)),(/nb,nx/))


end subroutine eigs


end module mod_calcss

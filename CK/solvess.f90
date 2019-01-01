! Main program for "Implications of Government Transfers for Aggregate Labor Market Fluctuations"
! by Youngsoo Jang, Takeki Sunakawa and Minchul Yum
! Updated on May 30, 2017 @ Pier's Cafe Todoroki
! Updated on June 8, 2017 @ Library for Social Sciences at Kobe U
! Updated on Jan 19, 2018 @ Mannheim: work decision based on reservation asset for each productivity
! Updated on Oct 17, 2018 @ Shanghai: some parameters are called via mod_parameters
! additional note: v4 eliminated tauk; found bug on Tax and transfers (fixed); introducing borrowing constraint
! Question: Y must be predicted? if we want the transfers to be a function of y/Yagg not y/Yss

program solvess


    use mod_utils
    use mod_parameters
    use mod_randnumber
    use mod_calcss
    use json_module
	implicit none


    real(8) agrid(na), mgrid(nm), xgrid(nx), pxij(nx,nx), mux(nx), ALOW
    real(8) bgrid(nb), pbij(nb,nb)
    integer idmat(ne,2), error
    real(8) w, r, zagg, Lagg, Kagg, Tragg, kgrid(nk), muk(nk), Ls, Ks, EPratio, wss, rss
    real(8) TrW5(5), EmpW5(5), stdlogwage, stdlogearn, rhoxest, Tr12, Tr1, Tr2, Tr45
    real(8) ndist(nk,ne), mumat(nk,ne), Wmat(na,ne), Nmat(na,ne), minasset(ne) !, mustationary(nk,ne) !, ndist2(nk,ne), mumat2(nk,ne), aa1, aa2, aa3, aa4
    integer ia, ik, im, ib, ix, jx, ie, iz !, iterout, iteroutvec(100), i, foreind, izvec(simTT), time, iv
    ! for json
    type(json_core) :: core
    type(json_value), pointer :: p, input, output, ss, irf
    ! integer, parameter :: ra = na-2
    ! integer, parameter :: rm = nm-2
    ! real(8) invTa(ra,ra), invTm(rm,rm)

    integer cto1, cto2, cro
    real(8) eptimeo


    call system_clock(cto1, cro)

    ! index for individual exogenous state
	do ib = 1,nb

		do ix = 1,nx

			ie = (ib-1)*nx + ix
			idmat(ie,1) = ib
			idmat(ie,2) = ix

		end do

	end do

    ! if (ssini0==1) then
    !     call readss(idmat,ALPHA,BETA,agrid,kgrid,bgrid,xgrid,pbij,pxij,xss,mumat,ndist,muk,Kss,Lss,Yss,ALOW)
    ! else
        write(*,*) "calculating the steady state..."
        zagg = 1.00d0
        print *, zagg
        call calcss(agrid,kgrid,bgrid,xgrid,pbij,pxij,mux,ALOW,idmat,1,0, &
    	        w,r,zagg,Lagg,Kagg,Tragg,mumat,ndist,Wmat,Nmat,minasset,Ls,Ks,EPratio,TrW5,EmpW5,stdlogwage,stdlogearn,rhoxest,error)
        ! call calcss(agrid,kgrid,bgrid,xgrid,pbij,pxij,xss,ALOW,idmat,1,1, &
  	    !     w,r,zagg,Lagg,Kagg,Tragg,mumat,ndist,Wmat,Nmat,minasset,Ls,Ks,EPratio,TrW5,EmpW5,stdlogwage,stdlogearn,rhoxest,error)
        ! pause
        !
        ! zagg = 1.00d0
        ! ! zagg = zagg+0.1d0*SDINOVZ
        ! ! zagg = zagg-3.0d0*SDINOVZ
        ! print *, zagg
        ! call calcss(agrid,kgrid,bgrid,xgrid,pbij,pxij,xss,ALOW,idmat,1,0, &
        ! 	w,r,zagg,Lagg,Kagg,Yagg,TrY,Trss,mumat2,ndist2,Ks,Ls,EPratio,TrW5,EmpW5,stdlogwage,stdlogearn,rhoxest,error)

        ! ! decomposition
        ! aa1 = 0.0d0
        ! aa2 = 0.0d0
        ! aa3 = 0.0d0
        ! aa4 = 0.0d0
        ! do ix = 1,nx
        !
        !     do ik = 1,nk
        !
        !         aa1 = aa1 + ndist1(ik,ix)*mumat1(ik,ix)
        !         aa2 = aa2 + ndist2(ik,ix)*mumat1(ik,ix)
        !         aa3 = aa3 + ndist1(ik,ix)*mumat2(ik,ix)
        !         aa4 = aa4 + ndist2(ik,ix)*mumat2(ik,ix)
        !
        !     end do
        !
        ! end do
        !
        ! print *, aa1, aa2, aa3, aa4

    ! end if
    ! write(*,"('  ALOW =', F10.5)") ALOW
    ! wss = (1.0d0-ALPHA)*Yss/Lss
    ! rss = ALPHA*Yss/Kss

    ! ! store the staionary distribution
    ! mustationary = mumat
    !
    ! ! set up the grid points
    ! mgrid = linspace((1.0d0-mpct)*Kss,(1.0d0+mpct)*Kss,nm)
    !
    ! ! grids for TFP
    ! call qnwnorm(nghz,0.0d0,SDINOVZ,xghz,wghz)
    ! if ((conzflag==2).or.(conzflag==3)) then
    !     ! zgrid = linspace(-mz*SDINOVZ/sqrt(1.0d0-RHOZ**2),mz*SDINOVZ/sqrt(1.0d0-RHOZ**2),nz)
    !     !call rouwenhorst(RHOZ,SDINOVZ,nz,zgrid,pzij)
    !     zgrid = linspace(-mz*SDINOVZ/sqrt(1.0d0-RHOZ**2),mz*SDINOVZ/sqrt(1.0d0-RHOZ**2),nz)
    ! else
    !     if (tauflagagg==1) then
    !         call tauchen(nz,0.0d0,RHOZ,SDINOVZ,mz,zgrid,pzij)
    !     else
    !         call rouwenhorst(RHOZ,SDINOVZ,nz,zgrid,pzij)
    !     end if
    ! end if
    !
    ! zgrid = exp(zgrid)
    !
    ! open(510, file="mgrid.txt")
    ! do im = 1,nm
    !     write(510,*) mgrid(im)
    ! end do
    ! close(510)
    !
    ! open(520, file="zgrid.txt")
    ! do iz = 1,nz
    !     write(520,*) zgrid(iz)
    ! end do
    ! close(520)
    !
    ! ! setup splines
    ! invTa = spbas(ra,agrid)
    ! invTm = spbas(rm,mgrid)
    ! ! take the inverse of the basis function matrix
    ! invTa = inv(invTa)
    ! invTm = inv(invTm)
    !
    ! !***************************************!
    ! !******* KRUSELL-SMITH ALGORITHM *******!
    ! !***************************************!
    ! ! generate random numbers
    ! izvec = genizvec(simTT,pzij)
    ! zvec = genzvec(simTT,RHOZ,SDINOVZ)
    !
    ! ! initial value function and forecasting rules
    ! Vmat = 0.0d0
    !
    ! if (conzflag==1) then
    !
    !     allocate(kappakp(2,nz),kappaw(2,nz),kappar(2,nz),kappal(2,nz),kappakpnew(2,nz),kappawnew(2,nz),kapparnew(2,nz),kappalnew(2,nz),rsquared(4,nz))
    !
    !     do iz = 1,nz
    !         ! NOTE: 031818 indices in kappas are changed (iz is now in columns)
    !         ! NOTE: The initial forecasting rules are from log-linearized representative-agent model
    !         ! First run with fcsteqn = 1 and outermkt = 0, then use the converged forecasting rule for r (and mp and w as well) by fcsteqn = 3
    !         !kappakp(1,iz) = (1.0d0-0.9431d0)*log(Kss) + 0.1359d0*log(zgrid(iz))
    !         !kappakp(2,iz) = 0.9431d0
    !         !kappaw(1,iz) = log(wss)-0.5714d0*log(Kss) + 0.3699d0*log(zgrid(iz))
    !         !kappaw(2,iz) = 0.5714d0
	! 		! NOTE: 23 May 2018: Following guesses work better
    !         kappakp(1,iz) = 0.12d0*(1.0d0+(zgrid(iz)-1.0d0)*0.4d0)
    !         kappakp(2,iz) = 0.95d0*(1.0d0+(zgrid(iz)-1.0d0)*0.01d0)
    !         kappaw(1,iz) = -0.23d0 + (zgrid(iz)-1.0d0)
    !         kappaw(2,iz) = 0.45d0/(1.0d0+(zgrid(iz)-1.0d0)*0.4d0)
    !         ! below is not used
    !         kappar(1,iz) = log(rss+DELTA)+1.0159d0*log(Kss) + 1.9813d0*log(zgrid(iz))
    !         kappar(2,iz) = -1.0159d0
    !         ! NOTE: How do we get this?
    !         ! based on baseline (near) solution
    !         kappal(1,iz) = -0.31d0 - (log(zgrid(iz)) - 1.0d0)
    !         kappal(2,iz) = -0.38d0 + 2.0d0*(log(zgrid(iz)) - 1.0d0)
    !
    !     end do
    !
    ! elseif (conzflag==3) then
    !
    !     allocate(kappakp(4,1),kappaw(4,1),kappar(4,1),kappal(4,1),kappakpnew(4,1),kappawnew(4,1),kapparnew(4,1),kappalnew(4,1),rsquared(4,1))
    !
    !     kappakp(1,1) = 0.113d0
    !     kappakp(2,1) = 0.954d0
    !     kappakp(3,1) = 0.100d0
    !     kappakp(4,1) = 0.0d0
    !     kappaw(1,1) = -0.398d0
    !     kappaw(2,1) = 0.522d0
    !     kappaw(3,1) = 0.800d0
    !     kappaw(4,1) = 0.0d0
    !     kappar(1,1) = -1.394d0
    !     kappar(2,1) = -0.799d0
    !     kappar(3,1) = 1.356d0
    !     kappar(4,1) = 0.0d0
    !     kappal(1,1) = -0.31d0
    !     kappal(2,1) = -0.799d0
    !     kappal(3,1) = 1.356d0
    !     kappal(4,1) = 0.0d0
    !
    ! else
    !
    !     allocate(kappakp(3,1),kappaw(3,1),kappar(3,1),kappal(3,1),kappakpnew(3,1),kappawnew(3,1),kapparnew(3,1),kappalnew(3,1),rsquared(4,1))
    !
    !     ! NOTE: how do we get these numbers?
    !     kappakp(1,1) = 0.113d0
    !     kappakp(2,1) = 0.954d0
    !     kappakp(3,1) = 0.100d0
    !     kappaw(1,1) = -0.398d0
    !     kappaw(2,1) = 0.522d0
    !     kappaw(3,1) = 0.800d0
    !     kappar(1,1) = -1.394d0
    !     kappar(2,1) = -0.799d0
    !     kappar(3,1) = 1.356d0
    !     kappal(1,1) = -0.31d0
    !     kappal(2,1) = -0.799d0
    !     kappal(3,1) = 1.356d0
    !
    ! end if
    !
    ! ! loading previous results
    ! if (vmatini0==1) call readwriteVmat(Vmat,1)
    ! if (fcstini0==1) call readwritekappa(kappakp,kappaw,kappar,kappal,ckappakp,ckappaw,ckappar,ckappal,conzflag,1) ! read
    !
    !
    ! !******* MAIN LOOP *******!
    ! diffout = 1d+4
    ! iterout = 1
    ! damp = dampout0
    ! bspct = bspct0
    !
    ! if (outermkt==1) foreind = 1  ! =1 flag for market clearing
    !
	! ! 25 May 2018: do while is replaced with do and if exit
    ! !do while (diffout>tolout)
	! do
    ! ! do iterout=1,3
    !
    !     ! NOTE: 030418 BE CAREFUL when use integer as a flag in if clause.  Need to specify value and use parenthesis (with .and., .or.)
    !
    !     !if ((iterout>1) .or. (vmatini0==0)) then
	! 	!if (vmatini0==0) then
    !         call inner(agrid,mgrid,invTa,invTm,bgrid,xgrid,zgrid,pbij,pxij,pzij,idmat,Yss,kappakp,kappaw,kappar,xghz,wghz,RHOZ,conzflag, &
    !             ALPHA,DELTA,HBAR,ALOW,B0,BETA,taul,wgtT,a0,a1,a2,Trscale,Transprog,T0, &
    !             tol,maxiter,diagnum,logdiff,vmatini,iterout,fcsteqn,linflag,howflag,Vmat,Wmat,Nmat,apmat,awpmat,anpmat,lmat)
    !
    !         call readwriteVmat(Vmat,0) ! write
    !     !end if
    !
    !     call outer(agrid,kgrid,mgrid,invTa,invTm,bgrid,xgrid,zgrid,pbij,pxij,pzij,simTT,izvec,zvec,idmat,Yss, &
    !         kappakp,kappaw,kappar,xghz,wghz,RHOZ,conzflag,Vmat,ALPHA,DELTA,HBAR,ALOW,B0,BETA,taul,wgtT,a0,a1,a2,Trscale,Transprog,T0, &
    !         tol,tolmkt,outermkt,foreind,bsctr,bspct,maxcountmkt,calcma,bsfixr,diagnumout,fcsteqn,linflag,ymat,mumat)
    !
    !     call calcforecast(drop,izvec,zgrid,ymat,DELTA,kappakpnew,kappawnew,kapparnew,kappalnew,rsquared,conzflag)
    !
    !     diffkp = maxval(maxval(abs(kappakpnew-kappakp),1),1)
    !     diffw = maxval(maxval(abs(kappawnew-kappaw),1),1)
    !     diffr = maxval(maxval(abs(kapparnew-kappar),1),1)
    !     diffl = maxval(maxval(abs(kappalnew-kappal),1),1)
    !
    !     if (fcsteqn==1) then
    !         diffout = max(diffkp,diffw)
    !     elseif (fcsteqn==2) then
    !         diffout = max(diffkp,diffr)
    !     else
    !         diffout = max(max(diffkp,diffw),diffr)
    !     end if
    !
    !     ! Once the forecasting rules are close to the ones in equiilibrium, the market-clearing conditions are strictly imposed.
    !     if ((outermkt==1) .and. (diffr<tolout*100.0d0)) then
    !         foreind = 1
    !     end if
    !
    !     iteroutvec(iterout) = iterout
    !     diffkpvec(iterout)  = diffkp
    !     diffwvec(iterout)   = diffw
    !     diffrvec(iterout)   = diffr
    !     difflvec(iterout)   = diffl
    !
    !     write(*,*) " "
    !     do i = 1,iterout
    !         write(*,"('  iteration ', I4, '   ||Tkappa-kappa|| = (', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ')')") &
    !         iteroutvec(i), diffkpvec(i), diffwvec(i), diffrvec(i), difflvec(i)
    !     end do
    !
    !     if (conzflag==1) then
    !
    !         do iz = 1,nz
    !             write(*,"('  z = ', F8.5, ': R2 of K = ', F8.5, ', R2 of w = ', F8.5, ', R2 of r = ', F8.5, ', R2 of L = ', F8.5)") &
    !             zgrid(iz), rsquared(1,iz), rsquared(2,iz), rsquared(3,iz), rsquared(4,iz)
    !         end do
    !         write(*,*) " "
    !
    !     else
    !
    !         write(*,"('  R2 of K = ', F8.5, ', R2 of w = ', F8.5, ', R2 of r = ', F8.5, ', R2 of L = ', F8.5)") &
    !             rsquared(1,1), rsquared(2,1), rsquared(3,1), rsquared(4,1)
    !         write(*,*) " "
    !
    !     end if
    !
    !     ! store the simulated time series and forecasting rules
    !     open(300, file="izvec.txt")
    !     open(310, file="mvec.txt")
    !     open(320, file="wvec.txt")
    !     open(330, file="rvec.txt")
    !     open(340, file="yvec.txt")
    !     open(350, file="cvec.txt")
    !     open(360, file="xvec.txt")
    !     open(370, file="lvec.txt")
    !     open(380, file="hvec.txt")
    !     open(390, file="zvec.txt")
    !     open(400, file="trvec.txt")
    !
    !     do time = 1,simT
    !
    !         write(300,*) izvec(time+drop)
    !         write(310,*) ymat(time+drop,1)
    !         write(320,*) ymat(time+drop,2)
    !         write(330,*) ymat(time+drop,3)
    !         write(340,*) ymat(time+drop,4)
    !         write(350,*) ymat(time+drop,5)
    !         write(360,*) ymat(time+drop,6)
    !         write(370,*) ymat(time+drop,7)
    !         write(380,*) ymat(time+drop,8)
    !         write(390,*) ymat(time+drop,9)
    !         write(400,*) ymat(time+drop,10)
    !
    !     end do
    !
    !     close(300)
    !     close(310)
    !     close(320)
    !     close(330)
    !     close(340)
    !     close(350)
    !     close(360)
    !     close(370)
    !     close(380)
    !     close(390)
    !     close(400)
    !
	! 	if (.NOT. diffout>tolout) exit
    !
    !     kappakp = (1.0d0-damp)*kappakpnew + damp*kappakp
    !     kappaw  = (1.0d0-damp)*kappawnew  + damp*kappaw
    !     kappar  = (1.0d0-damp)*kapparnew  + damp*kappar
    !     kappal  = (1.0d0-damp)*kappalnew  + damp*kappal
    !
	! 	! 24 May 2018 Note : saved forecasting rule should not be the updated one (at the final convergence)
	! 	call readwritekappa(kappakp,kappaw,kappar,kappal,ckappakp,ckappaw,ckappar,ckappal,conzflag,0) ! write
    !
	! 	iterout = iterout + 1
    !
    ! end do
    !
    !
    ! call system_clock(cto2, cro)
    ! eptimeo = (dble(cto2-cto1)/cro)/60.0d0
    ! write(*,"('  Total Elasped time (minutes) = ', F10.2)") eptimeo

    ! the i/o of data will be replaced by json
    ! output via json
    if (jsonoutput) then

        call core%initialize()
        call core%create_object(p, '')
        call core%create_object(input, 'input')
        call core%add(p, input)
        ! call core%add(input, 'nraflag', nraflag )
        call core%add(input, 'spblas', spblas )
        call core%add(input, 'linflag', linflag )
        call core%add(input, 'howflag', howflag )
        call core%add(input, 'logdiff', logdiff )
        ! call core%add(input, 'spliflag', spliflag )
        call core%add(input, 'transmat', transmat )
        ! call core%add(input, 'bisectm', bisectm )
        ! call core%add(input, 'fcstini', fcstini )
        ! call core%add(input, 'adjbias', adjbias )
        ! call core%add(input, 'naiveflag', naiveflag )
        call core%add(input, 'congrid', congrid )
        call core%add(input, 'deggrid', deggrid )
        ! call core%add(input, 'damp', damp )
        call core%add(input, 'dampss', dampss )
        call core%add(input, 'drop', drop )
        call core%add(input, 'simT', simT )
        ! call core%add(input, 'irfdrop', irfdrop )
        ! call core%add(input, 'irfT', irfT )
        call core%add(input, 'param', [ALPHA,DELTA,HBAR,RHOZ,SDINOVZ,taul,B0,BETA,RHOX,SDINOVX,T0,Trscale,Transprog,phi,ALOW] )
        call core%add(input, 'tol', [tolout,tol,toldist] )
        call core%add(input, 'agrid', agrid )
        call core%add(input, 'kgrid', kgrid )
        call core%add(input, 'bgrid', bgrid )
        call core%add(input, 'xgrid', xgrid )
        call core%add(input, 'pbij', reshape(pbij,(/nb*nb/)) )
        call core%add(input, 'pxij', reshape(pxij,(/nx*nx/)) )
        call core%add(input, 'mux', mux )
        nullify(input)

        call core%create_object(ss, 'ss')
        call core%add(p, ss)
        call core%add(ss, 'Wmat', reshape(Wmat,(/na*ne/)) )
        call core%add(ss, 'Nmat', reshape(Nmat,(/na*ne/)) )
        call core%add(ss, 'mumat', reshape(mumat,(/nk*ne/)) )
        call core%add(ss, 'ndist', reshape(ndist,(/nk*ne/)) )
        call core%add(ss, 'minasset', minasset )
        ! call core%add(ss, 'psix',   reshape(psix,(/nx*nz/)) )
        ! call core%add(ss, 'mnow',   mnow )
        ! call core%add(ss, 'shareWvec',     ShareWvec )
        ! call core%add(ss, 'gini',   gini )
        ! call core%add(ss, 'thresholdWvec', ThresholdWvec )
        nullify(ss)
        call core%print(p, trim(jsonfilename))
        call core%destroy(p)
        print *, trim(jsonfilename)

        ! call core%create_object(output, 'output')
        ! call core%add(p, output)
        ! call core%add(output, 'vmat0',  reshape(vmat0,(/nk*nm*nz*nx/)) )
        ! call core%add(output, 'gmat0',  reshape(gmat0,(/nk*nm*nz*nx/)) )
        ! call core%add(output, 'mpmat0', reshape(mpmat0,(/nm*nz/)) )
        !
        ! call core%add(output, 'Yvec',   aggsim(:,1) )
        ! call core%add(output, 'Zvec',   aggsim(:,2) )
        ! call core%add(output, 'Nvec',   aggsim(:,3) )
        ! call core%add(output, 'Cvec',   aggsim(:,4) )
        ! call core%add(output, 'Ivec',   aggsim(:,5) )
        ! call core%add(output, 'Kvec',   aggsim(:,6) )
        ! call core%add(output, 'Kpvec',  aggsim(:,7) )
        ! call core%add(output, 'Xvec',   aggsim(:,8) )
        ! call core%add(output, 'shareW1',    disaggsim(:,1) )
        ! call core%add(output, 'shareW2',    disaggsim(:,2) )
        ! call core%add(output, 'shareW3',    disaggsim(:,3) )
        ! call core%add(output, 'shareW4',    disaggsim(:,4) )
        ! call core%add(output, 'shareW5',    disaggsim(:,5) )
        ! call core%add(output, 'shareW9095', disaggsim(:,6) )
        ! call core%add(output, 'shareW9599', disaggsim(:,7) )
        ! call core%add(output, 'shareWT1',   disaggsim(:,8) )
        ! call core%add(output, 'gini',       disaggsim(:,9) )
        !
        ! call core%add(output, 'eptimein',   epvec(1:iter,1) )
        ! call core%add(output, 'eptimeout',  epvec(1:iter,2) )
        ! nullify(output)
        !
        ! call core%create_object(irf, 'irf')
        ! call core%add(p, irf)
        ! call core%add(irf, 'Yvec', aggirf(:,1) )
        ! call core%add(irf, 'Zvec', aggirf(:,2) )
        ! call core%add(irf, 'Nvec', aggirf(:,3) )
        ! call core%add(irf, 'Cvec', aggirf(:,4) )
        ! call core%add(irf, 'Ivec', aggirf(:,5) )
        ! call core%add(irf, 'Kvec', aggirf(:,6) )
        ! call core%add(irf, 'shareW1',    disaggirf(:,1) )
        ! call core%add(irf, 'shareW2',    disaggirf(:,2) )
        ! call core%add(irf, 'shareW3',    disaggirf(:,3) )
        ! call core%add(irf, 'shareW4',    disaggirf(:,4) )
        ! call core%add(irf, 'shareW5',    disaggirf(:,5) )
        ! call core%add(irf, 'shareW9095', disaggirf(:,6) )
        ! call core%add(irf, 'shareW9599', disaggirf(:,7) )
        ! call core%add(irf, 'shareWT1',   disaggirf(:,8) )
        ! call core%add(irf, 'gini',       disaggirf(:,9) )
        ! call core%add(irf, 'mu0',  reshape(mu0,(/nb*nx/)) )
        ! nullify(irf)
        !
        ! call core%print(p, trim(jsonfilename))
        ! call core%destroy(p)
        ! print *, trim(jsonfilename)

    end if


end program solvess

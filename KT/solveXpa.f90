program solveXpa


    use mod_parameters
    use mod_functions
    use mod_spline
    use mod_calcss
    use mod_planner
    use mod_inner
    use mod_outer
    use mod_simulation
    use json_module
	implicit none


    real(8) kbounds(2), mbounds(2)
    real(8) Gz(nz), Pz(nz,nz), Ge(ne), Pe(ne,ne), mue(ne), knotsk(nk), knotsb(nb), invTk(rk,rk), knotsm(nm), invTm(rm,rm)
    real(8) muss(nb,ne), mu0(nb,ne), evec(ne,2), mnow, psie(ne), mpmat(nb,ne), ymat(nb,ne), nmat(nb,ne), ikvec(7), mp, ynow, nnow, inow
    integer iz, jz, ie, je, ik, im, ib, iter, ct1, ct2, cr, tt, error
    real(8) vmatss(nk,ne), gmatss(nk,ne), nmatpl(nk,nz), vmatpl(nk,nz), gmatpl(nk,nz), mpmat0(nm,nz), pmat0(nm,nz), mpmat1(nm,nz), pmat1(nm,nz)
    real(8) vmat0(nk,nm,nz,ne), gmat0(nk,nm,nz,ne)
    real(8) diff, diffmp, diffp, eptimein, eptimeout, eptime, epvec(100,2)
    real(8) cumPz(nz,nz), rand, aggsim(simT+drop,8), disaggsim(simT+drop,7)
    real(8) Zvec(irfTT), aggirf(irfTT,8), disaggirf(irfTT,7)
    integer izvec(simT+drop)
    integer seedsize
    integer, allocatable :: seed(:)
    ! for json
    type(json_file) :: json
    type(json_core) :: core
    real(8), allocatable :: temp(:)
    integer found
    type(json_value), pointer :: p, input, output, ss, irf
    ! for dgetrf, dgetri
    integer INFO
    integer IPIVk(rk), IPIVm(rm)
    real(8) WORKk(rk), WORKm(rm)


    call system_clock(ct1, cr)

    if (nz==1) then
        Gz = 1.0d0
        Pz = 1.0d0
    else
        call tauchen(nz,0.0d0,RHO,SIGMA,mz,Gz,Pz)
        Gz = exp(Gz)
    end if

    if (ne==1) then
        Ge = 1.0d0
        Pe = 1.0d0
        mue = 1.0d0
    else
        call tauchen(ne,0.0d0,RHOE,SIGE,me,Ge,Pe)
        Ge = exp(Ge)
        mue = calcmu(Pe)
    end if

    kbounds = (/0.1d0, 5.0d0/)

    knotsk = logspace(log(kbounds(1)-kbounds(1)+1.0d0)/log(10.0d0), log(kbounds(2)-kbounds(1)+1.0d0)/log(10.0d0), nk)
    knotsk = knotsk + (kbounds(1) - 1.0d0)
    knotsb = linspace(kbounds(1), kbounds(2), nb)
    invTk = spbas(rk,knotsk)
    call dgetrf(rk,rk,invTk,rk,IPIVk,INFO)
    call dgetri(rk,invTk,rk,IPIVk,WORKk,rk,INFO)

    print *, 'Calculating the steady state'
    call calcss(knotsk,knotsb,invTk,Ge,Pe,mue,vmatss,gmatss,muss,evec,mpmat,ymat,nmat)
    mu0 = muss
    call calcdiststat(knotsb,mu0,mpmat,ikvec)
    print *, ikvec
    if (adjbias .eqv. .false.) evec = 0.0d0

    ! fraction of e-indexed capital
    mnow = 0.0d0
    do ie = 1,ne
        psie(ie) = sum(knotsb*mu0(:,ie),1)/sum(mu0(:,ie),1)
        mnow = mnow + sum(knotsb*mu0(:,ie),1)
    end do
    psie = psie/mnow
    print *, evec, psie, mnow

    mbounds = (/0.75d0*mnow, 1.25d0*mnow/)

    knotsm = linspace(mbounds(1), mbounds(2), nm)
    invTm = spbas(rm,knotsm)
    call dgetrf(rm,rm,invTm,rm,IPIVm,INFO)
    call dgetri(rm,invTm,rm,IPIVm,WORKm,rm,INFO)

    print *, 'Solving the planner''s problem'
    call solvepl(knotsk,invTk,Gz,Pz,gmatpl,nmatpl,vmatpl)
    call forecastpl(nmatpl,knotsk,knotsm,Gz,mpmat0,pmat0)


    diff = 1d+4
    iter = 0

    print *, 'Solving the lumpy model'
    do while (diff>critout)

        call inner(mpmat0,pmat0,knotsk,knotsm,invTk,invTm,Gz,Pz,Ge,Pe,vmat0,gmat0,eptimein,error)

        call outer(vmat0,mpmat0,knotsk,knotsm,invTk,invTm,Gz,Pz,Ge,Pe,psie,mue,evec,mpmat1,pmat1,eptimeout)

        diffmp = maxval(maxval(abs(mpmat1-mpmat0),1),1)
        diffp = maxval(maxval(abs(pmat1-pmat0),1),1)
        diff = max(diffmp,diffp)
        iter = iter + 1
        ! diagnosis
        write(*,"('  iteration ', I4, '  ||Tmp-mp|| = ', F10.5, '  ||Tp-p|| = ', F10.5)") iter, diffmp, diffp

        mpmat0 = damp*mpmat1 + (1.0d0-damp)*mpmat0
        pmat0  = damp*pmat1 + (1.0d0-damp)*pmat0

        epvec(iter,1) = eptimein
        epvec(iter,2) = eptimeout

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime
    print *, sum(epvec(1:iter,:),1)/iter


    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed)
    ! print *, "Size of seed array is", seedsize
    call random_seed(put=seed)

    cumPz = cumsum(Pz)
    izvec(1) = ceiling(dble(nz)/2.0d0)
    do tt = 1,simT+drop-1

        call random_number(rand)
        izvec(tt+1) = count(rand-cumPz(izvec(tt),:)>=0)
        izvec(tt+1) = min(izvec(tt+1)+1,nz)

    end do

    print *, 'Simulating the lumpy model'
    call simulation(vmat0,mpmat0,pmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Pz,Ge,Pe,izvec,aggsim,disaggsim,mu0,eptimeout)


    ! determine the path of TFP shock
    Zvec(1:irfdrop+1) = 1.0d0
    Zvec(irfdrop+2) = max(exp(RHO*log(Zvec(irfdrop+1))+shocksize),Gz(1))

	do tt=irfdrop+2,irfTT-1
        Zvec(tt+1) = exp(RHO*log(Zvec(tt)))
    end do

    call calcirf1(vmat0,mpmat0,pmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Pz,Ge,Pe,Zvec,aggirf,disaggirf,mu0,eptimeout)


    ! output via json
    if (jsonoutput) then

        call core%initialize()
        call core%create_object(p, '')
        call core%create_object(input, 'input')
        call core%add(p, input)
        call core%add(input, 'nraflag', nraflag )
        call core%add(input, 'bisectw', bisectw )
        call core%add(input, 'outermkt', outermkt )
        call core%add(input, 'fcstini', fcstini )
        call core%add(input, 'adjbias', adjbias )
        call core%add(input, 'naiveflag', naiveflag )
        call core%add(input, 'damp', damp )
        call core%add(input, 'dampss', dampss )
        call core%add(input, 'drop', drop )
        call core%add(input, 'simT', simT )
        call core%add(input, 'irfdrop', irfdrop )
        call core%add(input, 'irfT', irfT )
        call core%add(input, 'param', [GAMY,BETA,DELTA,NU,RHO,SIGMA,THETA,ETA,RHOE,SIGE] )
        call core%add(input, 'crit', [critout,critmu,critin,critbp,critbn,critg,critn] )
        call core%add(input, 'knotsk', knotsk )
        call core%add(input, 'knotsm', knotsm )
        call core%add(input, 'knotsb', knotsb )
        call core%add(input, 'Gz', Gz )
        call core%add(input, 'Pz', reshape(Pz,(/nz*nz/)) )
        call core%add(input, 'Ge', Ge )
        call core%add(input, 'Pe', reshape(Pe,(/ne*ne/)) )
        call core%add(input, 'mue', mue )
        call core%add(input, 'izvec', izvec )
        nullify(input)

        call core%create_object(ss, 'ss')
        call core%add(p, ss)
        call core%add(ss, 'vmatss', reshape(vmatss,(/nk*ne/)) )
        call core%add(ss, 'gmatss', reshape(gmatss,(/nk*ne/)) )
        call core%add(ss, 'muss',   reshape(muss,(/nb*ne/)) )
        call core%add(ss, 'mpmat',  reshape(mpmat,(/nb*ne/)) )
        call core%add(ss, 'ymat',   reshape(ymat,(/nb*ne/)) )
        call core%add(ss, 'nmat',   reshape(nmat,(/nb*ne/)) )
        call core%add(ss, 'evec',   reshape(evec,(/ne*2/)) )
        call core%add(ss, 'psie',   psie )
        call core%add(ss, 'mnow',   mnow )
        nullify(ss)
        ! call core%print(p, trim(jsonfilename))
        ! call core%destroy(p)
        ! print *, trim(jsonfilename)
        ! pause

        call core%create_object(output, 'output')
        call core%add(p, output)
        call core%add(output, 'vmat0',  reshape(vmat0,(/nk*nm*nz*ne/)) )
        call core%add(output, 'gmat0',  reshape(gmat0,(/nk*nm*nz*ne/)) )
        call core%add(output, 'mpmat0', reshape(mpmat0,(/nm*nz/)) )
        call core%add(output, 'pmat0',  reshape(pmat0,(/nm*nz/)) )
        call core%add(output, 'Yvec',   aggsim(:,1) )
        call core%add(output, 'Zvec',   aggsim(:,2) )
        call core%add(output, 'Nvec',   aggsim(:,3) )
        call core%add(output, 'Cvec',   aggsim(:,4) )
        call core%add(output, 'Ivec',   aggsim(:,5) )
        call core%add(output, 'Kvec',   aggsim(:,6) )
        call core%add(output, 'Kpvec',  aggsim(:,7) )
        call core%add(output, 'Xvec',   aggsim(:,8) )
        call core%add(output, 'ikmean',     disaggsim(:,1) )
        call core%add(output, 'ikstddev',   disaggsim(:,2) )
        call core%add(output, 'ikinaction', disaggsim(:,3) )
        call core%add(output, 'ikspikepos', disaggsim(:,4) )
        call core%add(output, 'ikspikeneg', disaggsim(:,5) )
        call core%add(output, 'ikpos',      disaggsim(:,6) )
        call core%add(output, 'ikneg',      disaggsim(:,7) )
        call core%add(output, 'eptimein',   epvec(1:iter,1) )
        call core%add(output, 'eptimeout',  epvec(1:iter,2) )
        nullify(output)

        call core%create_object(irf, 'irf')
        call core%add(p, irf)
        call core%add(irf, 'Yvec', aggirf(:,1) )
        call core%add(irf, 'Zvec', aggirf(:,2) )
        call core%add(irf, 'Nvec', aggirf(:,3) )
        call core%add(irf, 'Cvec', aggirf(:,4) )
        call core%add(irf, 'Ivec', aggirf(:,5) )
        call core%add(irf, 'Kvec', aggirf(:,6) )
        call core%add(irf, 'ikmean',     disaggirf(:,1) )
        call core%add(irf, 'ikstddev',   disaggirf(:,2) )
        call core%add(irf, 'ikinaction', disaggirf(:,3) )
        call core%add(irf, 'ikspikepos', disaggirf(:,4) )
        call core%add(irf, 'ikspikeneg', disaggirf(:,5) )
        call core%add(irf, 'ikpos',      disaggirf(:,6) )
        call core%add(irf, 'ikneg',      disaggirf(:,7) )
        nullify(irf)

        call core%print(p, trim(jsonfilename))
        call core%destroy(p)
        print *, trim(jsonfilename)

    end if


end program solveXpa

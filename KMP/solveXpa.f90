program solveXpa


    use mod_functions
    use mod_parameters
    use mod_spline
    use mod_calcss
    use mod_inner
    use mod_outer
    use json_module
    implicit none

    real(8) kbounds(2), Gz(nz), Gl(nz), Pz(nz,nz), Gez(ne,nz), Pez(ne,ne,2*nz), knotsk(nk), knotsm(nm), knotsb(nb), invTk(rk,rk), invTm(rm,rm)
    real(8) Gy(ny), Py(ny,ny), Gd(nd), Pd(nd,nd)
    integer ik, ib, ix, ie, iy, id, im, iz, idmat(nx,3), iter, ct1, ct2, cr, tt
    real(8) vmatss(nk,nx), gmatss(nk,nx), muss(nb,nx), mu0(nb,nx), evec(nx,nz), mnow, psix(nx,nz), mpmat(nb,nx)
    real(8) znow, lnow, weightGE, weightBE, weightGU, weightBU, Ge(ne), Pe(ne,ne), muez(ne,nz), muy(ny), muxz(nx,nz), mux(nx)
    integer ThresholdWvec(7)
    real(8) ShareWvec(8), gini
    real(8) BetaK(nz,2), R2(nz), mpmat0(nm,nz), mpmat1(nm,nz), vmat0(nk,nm,nz,nx), gmat0(nk,nm,nz,nx)
    real(8) diff, eptime, eptimein, eptimeout, epvec(100,2)
    real(8) cumPz(nz,nz), rand, aggsim(simT+drop,8), disaggsim(simT+drop,9), aggirf(irfT+irfdrop,8), disaggirf(irfT+irfdrop,9)
    integer izvec(simT+drop), izvecirf(irfT+irfdrop)
    integer seedsize
    integer, allocatable :: seed(:)
    ! for json
    type(json_core) :: core
    type(json_value), pointer :: p, input, output, ss, irf
    ! for dgetrf, dgetri
    integer INFOk, INFOm
    integer IPIVk(rk), IPIVm(rm)
    real(8) WORKk(rk), WORKm(rm)


    call system_clock(ct1,cr)

    ! endogenous variables' grid points
    kbounds = (/0.0d0, kmax/)

    if (congrid==1) then
        ! more grid points concentrated toward the left
        knotsk = logspace(0.0d0,deggrid,nk) - 1.0d0
        knotsk = knotsk/(10.0d0**deggrid - 1.0d0)
        knotsk = knotsk*(kbounds(2)-kbounds(1)) + kbounds(1)
        knotsb = logspace(0.0d0,deggrid,nb) - 1.0d0
        knotsb = knotsb/(10.0d0**deggrid - 1.0d0)
        knotsb = knotsb*(kbounds(2)-kbounds(1)) + kbounds(1)
    else
        knotsk = logspace(log(kbounds(1) - 1.0d0*kbounds(1) + 1.0d0)/log(10.0d0), log(kbounds(2) - 1.0d0*kbounds(1) + 1.0d0)/log(10.0d0), nk)
        knotsk = knotsk + kbounds(1) - 1.0d0
        ! knotsb = linspace(kbounds(1),kbounds(2),nb)
        knotsb = logspace(log(kbounds(1) - 1.0d0*kbounds(1) + 1.0d0)/log(10.0d0), log(kbounds(2) - 1.0d0*kbounds(1) + 1.0d0)/log(10.0d0), nb)
        knotsb = knotsb + kbounds(1) - 1.0d0
    end if

    invTk = spbas(rk,knotsk)
    call dgetrf(rk,rk,invTk,rk,IPIVk,INFOk)
    call dgetri(rk,invTk,rk,IPIVk,WORKk,rk,INFOk)

    ! ! Aiyagari
    ! Ge(1) = mu
    ! Ge(2) = 1.0d0
    ! Pe(1,1) = puu
    ! Pe(1,2) = 1.0d0-puu
    ! Pe(2,1) = 1.0d0-pee
    ! Pe(2,2) = pee
    ! Gy = 1.0d0
    ! Py = 1.0d0
    ! Gd = 0.96d0
    ! Pd = 1.0d0
    ! ! KS 1998: no tax for unemployed
    ! Gez(1,1) = mu
    ! Gez(2,1) = h*(1.0d0-taug)
    ! Gez(1,2) = mu
    ! Gez(2,2) = h*(1.0d0-taub)
    ! call calcPeKS(pgg,pbb,ug,ub,DurationUg,DurationUb,corr,Pez)
    ! Gy = 1.0d0
    ! Py = 1.0d0
    ! Gd = 0.989975d0
    ! Pd = 1.0d0
    ! KMP: tax for both
    Gez(1,1) = mu*(1.0d0-taug)
    Gez(2,1) = h*(1.0d0-taug)
    Gez(1,2) = mu*(1.0d0-taub)
    Gez(2,2) = h*(1.0d0-taub)
    call calcPeKMP(pgg,pbb,pgg00,pgg10,pgb00,pgb10,pbg00,pbg10,pbb00,pbb10,rhoy,sigy,dbar,deps,Pez,Gy,Py,Gd,Pd)

    ! index for exogenous grid
    do id = 1,nd

        do iy = 1,ny

            do ie = 1,ne

                ix = (id-1)*ny*ne+(iy-1)*ne+ie
                idmat(ix,1) = ie
                idmat(ix,2) = iy
                idmat(ix,3) = id

            end do

        end do

    end do

    ! exogenous grid for steady state
    ! NOTE: weightGU+weightBU=1 and weightGE+weightBE=1 hold
    weightGE = (1.0d0-FractionZb) * (1.0d0-ug)/(1.0d0-uss)
    weightBE = FractionZb         * (1.0d0-ub)/(1.0d0-uss)
    weightGU = (1.0d0-FractionZb) * ug/uss
    weightBU = FractionZb         * ub/uss
    Ge(1) = weightGU*Gez(1,1) + weightBU*Gez(1,2)
    Ge(2) = weightGE*Gez(2,1) + weightBE*Gez(2,2)
    Pe(1,1) = weightGU*(Pez(1,1,1) + Pez(1,1,2)) + weightBU*(Pez(1,1,3) + Pez(1,1,4))
    Pe(1,2) = weightGU*(Pez(1,2,1) + Pez(1,2,2)) + weightBU*(Pez(1,2,3) + Pez(1,2,4))
    Pe(2,1) = weightGE*(Pez(2,1,1) + Pez(2,1,2)) + weightBE*(Pez(2,1,3) + Pez(2,1,4))
    Pe(2,2) = weightGE*(Pez(2,2,1) + Pez(2,2,2)) + weightBE*(Pez(2,2,3) + Pez(2,2,4))
    ! Ge = (/0.150000000000000, 1.06397849462366/)
    ! Pe(1,:) = (/0.507440476190476, 0.492559523809524/)
    ! Pe(2,:) = (/3.707437275985664E-002, 0.962925627240143/)
    print *, Ge
    print *, Pe(1,:), sum(Pe(1,:))
    print *, Pe(2,:), sum(Pe(2,:))
    ! pause

    ! stationary distribution for exogenous grid
    muez(1,:) = (/ug, ub/)
    muez(2,:) = (/1.0d0-ug, 1.0d0-ub/)

    if (ny==1 .and. nd==1) then

        muxz(1:ne,:) = muez

    else

        muy = calcmu(Py)

        do id=1,nd

            do iy=1,ny

                muxz(ny*ne*(id-1)+ne*(iy-1)+1:ny*ne*(id-1)+ne*iy,1) = muez(:,1)*muy(iy)*(1.0d0/dble(nd))
                muxz(ny*ne*(id-1)+ne*(iy-1)+1:ny*ne*(id-1)+ne*iy,2) = muez(:,2)*muy(iy)*(1.0d0/dble(nd))

            end do

        end do

    end if

    ! mux = (1.0d0-FractionZb)*muxz(:,1) + FractionZb*muxz(:,2)
    mux = calcmu(Pe)
    print *, calcmu(Pez(:,:,1)/pgg)
    print *, calcmu(Pez(:,:,2)/(1.0d0-pgg))
    print *, calcmu(Pez(:,:,3)/(1.0d0-pbb))
    print *, calcmu(Pez(:,:,4)/pbb)
    print *, uss, 1.0d0-uss
    print *, mux
    ! pause

    print *, 'Calculating the steady state'
    znow = 1.0d0
    ! lnow = Ge(1)*mux(1)+Ge(2)*mux(2) !h*(1.0d0-uss)
    lnow = h*(1.0d0-uss)
    print *, lnow
    ! pause
    call calcss(knotsk,knotsb,invTk,znow,lnow,Ge,Gy,Gd,Pe,Py,Pd,idmat,muxz,vmatss,gmatss,muss,evec,mpmat)
    mu0 = muss
    call calcdiststat(knotsb,mu0,ThresholdWvec,ShareWvec,gini)
    ! print *, knotsb(ThresholdWvec)
    print *, ShareWvec, gini
    if (adjbias .eqv. .false.) evec = 0.0d0

    ! fraction of e-indexed capital
    mnow = 0.0d0
    do ix = 1,nx
        psix(ix,1) = sum(knotsb*mu0(:,ix),1)/muxz(ix,1)
        psix(ix,2) = sum(knotsb*mu0(:,ix),1)/muxz(ix,2)
        mnow = mnow + sum(knotsb*mu0(:,ix),1)
    end do
    psix = psix/mnow
    print *, mnow, psix, evec
    ! pause

    ! exogenous variables' grid points
    Gz = (/zg, zb/)
    Pz(1,:) = (/pgg, 1.0d0-pgg/)
    Pz(2,:) = (/1.0d0-pbb, pbb/)
    Gl = (/h*(1.0d0-ug), h*(1.0d0-ub)/)

    ! grid for aggregate capital
    knotsm = linspace(0.80d0*mnow, 1.15d0*mnow, nm)
    invTm = spbas(rm,knotsm)
    call dgetrf(rm,rm,invTm,rm,IPIVm,INFOm)
    call dgetri(rm,invTm,rm,IPIVm,WORKm,rm,INFOm)

    ! initial forecasting rules
    BetaK = 0.0d0
    ! from MMV2010JEDC
    BetaK(1,1) = 0.137800d0
    BetaK(1,2) = 0.963238d0
    BetaK(2,1) = 0.123815d0
    BetaK(2,2) = 0.965565d0
    call Beta2mat(BetaK,knotsm,mpmat0)


    diff = 1d+4
    iter = 0

    print *, 'Solving the model'
    do while(diff>critout)

        call inner(mpmat0,knotsk,knotsm,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,iter,vmat0,gmat0,eptimein)

        call outer(vmat0,gmat0,mpmat0,knotsk,knotsm,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,psix,muxz,evec,mpmat1,eptimeout)

        diff = maxval(maxval(abs(mpmat1-mpmat0),1),1)
        iter = iter + 1
        ! diagnosis
        write(*,"('  iteration ', I4, '  ||Tmp-mp|| = ', F10.5)") iter, diff

        mpmat0 = damp*mpmat1 + (1.0d0-damp)*mpmat0

        epvec(iter,1) = eptimein
        epvec(iter,2) = eptimeout

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime
    print *, sum(epvec(1:iter,:),1)/iter

    print *, knotsm(1), knotsm(nm)
    print *, minval(minval(mpmat0,1),1), maxval(maxval(mpmat0,1),1)

    ! simulation
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed)
    ! print *, "Size of seed array is", seedsize
    call random_seed(put=seed)

    cumPz = cumsum(Pz)
    izvec(1) = 1 !ceiling(dble(nz)/2.0d0)
    do tt = 1,simT+drop-1

        call random_number(rand)
        izvec(tt+1) = count(rand-cumPz(izvec(tt),:)>=0.0d0)
        izvec(tt+1) = min(izvec(tt+1)+1,nz)

    end do

    print *, 'Simulating the model'
    call simulation(vmat0,gmat0,mpmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,izvec,aggsim,disaggsim,mu0,eptimeout)

    print *, 'Simulating the model with impulse responses'
    izvecirf = 1
    izvecirf(irfdrop+2) = 2
    mu0 = muss
    call simulation(vmat0,gmat0,mpmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,izvecirf,aggirf,disaggirf,mu0,eptimeout)
    print *, 'done.'


    ! the i/o of data will be replaced by json
    ! output via json
    if (jsonoutput) then

        call core%initialize()
        call core%create_object(p, '')
        call core%create_object(input, 'input')
        call core%add(p, input)
        call core%add(input, 'nraflag', nraflag )
        call core%add(input, 'linflag', linflag )
        call core%add(input, 'spliflag', spliflag )
        call core%add(input, 'transmat', transmat )
        call core%add(input, 'bisectm', bisectm )
        ! call core%add(input, 'fcstini', fcstini )
        call core%add(input, 'adjbias', adjbias )
        call core%add(input, 'naiveflag', naiveflag )
        ! call core%add(input, 'congrid', congrid )
        ! call core%add(input, 'deggrid', deggrid )
        call core%add(input, 'damp', damp )
        call core%add(input, 'dampss', dampss )
        call core%add(input, 'drop', drop )
        call core%add(input, 'simT', simT )
        call core%add(input, 'irfdrop', irfdrop )
        call core%add(input, 'irfT', irfT )
        call core%add(input, 'param', [THETA,SIGMA,ALPHA,DELTA,zg,zb,ug,ub,mu,h,taug,taub,rhoy,sigy,dbar,deps] )
        call core%add(input, 'crit', [critout,critbp,critin,critg,critn,critmu] )
        call core%add(input, 'knotsk', knotsk )
        call core%add(input, 'knotsm', knotsm )
        call core%add(input, 'knotsb', knotsb )
        call core%add(input, 'Gz', Gz )
        call core%add(input, 'Gl', Gl )
        call core%add(input, 'Pz', reshape(Pz,(/nz*nz/)) )
        call core%add(input, 'Ge', Ge )
        call core%add(input, 'Pez', reshape(Pez,(/ne*ne*nz/)) )
        call core%add(input, 'Pe', reshape(Pe,(/ne*ne/)) )
        call core%add(input, 'Gy', Gy )
        call core%add(input, 'Py', reshape(Py,(/ny*ny/)) )
        call core%add(input, 'Gd', Gd )
        call core%add(input, 'Pd', reshape(Pd,(/nd*nd/)) )
        call core%add(input, 'muez', reshape(muez, (/ne*nz/)) )
        call core%add(input, 'muy', muy )
        call core%add(input, 'muxz', reshape(muxz, (/nx*nz/)) )
        call core%add(input, 'mux', mux )
        call core%add(input, 'izvec', izvec )
        nullify(input)

        call core%create_object(ss, 'ss')
        call core%add(p, ss)
        call core%add(ss, 'vmatss', reshape(vmatss,(/nk*ne/)) )
        call core%add(ss, 'gmatss', reshape(gmatss,(/nk*ne/)) )
        call core%add(ss, 'muss',   reshape(muss,(/nb*nx/)) )
        call core%add(ss, 'mpmat',  reshape(mpmat,(/nb*nx/)) )
        call core%add(ss, 'evec',   reshape(evec,(/nx*nz/)) )
        call core%add(ss, 'psix',   reshape(psix,(/nx*nz/)) )
        call core%add(ss, 'mnow',   mnow )
        call core%add(ss, 'shareWvec',     ShareWvec )
        call core%add(ss, 'gini',   gini )
        ! call core%add(ss, 'thresholdWvec', ThresholdWvec )
        nullify(ss)
        ! call core%print(p, trim(jsonfilename))
        ! call core%destroy(p)
        ! print *, trim(jsonfilename)

        call core%create_object(output, 'output')
        call core%add(p, output)
        call core%add(output, 'vmat0',  reshape(vmat0,(/nk*nm*nz*nx/)) )
        call core%add(output, 'gmat0',  reshape(gmat0,(/nk*nm*nz*nx/)) )
        call core%add(output, 'mpmat0', reshape(mpmat0,(/nm*nz/)) )

        call core%add(output, 'Yvec',   aggsim(:,1) )
        call core%add(output, 'Zvec',   aggsim(:,2) )
        call core%add(output, 'Nvec',   aggsim(:,3) )
        call core%add(output, 'Cvec',   aggsim(:,4) )
        call core%add(output, 'Ivec',   aggsim(:,5) )
        call core%add(output, 'Kvec',   aggsim(:,6) )
        call core%add(output, 'Kpvec',  aggsim(:,7) )
        call core%add(output, 'Xvec',   aggsim(:,8) )
        call core%add(output, 'shareW1',    disaggsim(:,1) )
        call core%add(output, 'shareW2',    disaggsim(:,2) )
        call core%add(output, 'shareW3',    disaggsim(:,3) )
        call core%add(output, 'shareW4',    disaggsim(:,4) )
        call core%add(output, 'shareW5',    disaggsim(:,5) )
        call core%add(output, 'shareW9095', disaggsim(:,6) )
        call core%add(output, 'shareW9599', disaggsim(:,7) )
        call core%add(output, 'shareWT1',   disaggsim(:,8) )
        call core%add(output, 'gini',       disaggsim(:,9) )

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
        call core%add(irf, 'shareW1',    disaggirf(:,1) )
        call core%add(irf, 'shareW2',    disaggirf(:,2) )
        call core%add(irf, 'shareW3',    disaggirf(:,3) )
        call core%add(irf, 'shareW4',    disaggirf(:,4) )
        call core%add(irf, 'shareW5',    disaggirf(:,5) )
        call core%add(irf, 'shareW9095', disaggirf(:,6) )
        call core%add(irf, 'shareW9599', disaggirf(:,7) )
        call core%add(irf, 'shareWT1',   disaggirf(:,8) )
        call core%add(irf, 'gini',       disaggirf(:,9) )
        call core%add(irf, 'mu0',  reshape(mu0,(/nb*nx/)) )
        nullify(irf)

        call core%print(p, trim(jsonfilename))
        call core%destroy(p)
        print *, trim(jsonfilename)

    end if


end program solveXpa

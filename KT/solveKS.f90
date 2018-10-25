program solveKS


    use mod_parameters
    use mod_functions
    use mod_spline
    use mod_calcss
    use mod_planner
    use mod_inner
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
    real(8) BetaK(nz,2), Betap(nz,2), R2(nz,2)
    real(8) cumPz(nz,nz), rand, aggsim(simT+drop,8), disaggsim(simT+drop,7)
    real(8) Kvec(simT+drop), Cvec(simT+drop) !, wm, cnow, mvec(simT+drop), pvec(simT+drop), DHmax(nz,2), DHmean(nz,2)
    integer izvec(simT+drop)
    integer seedsize
    integer, allocatable :: seed(:)
    ! for json
    type(json_file) :: json
    type(json_core) :: core
    real(8), allocatable :: temp(:)
    integer found
    type(json_value), pointer :: p, input, output, ss
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
    ! kbounds = (/0.1d0, 8.0d0/)

    knotsk = logspace(log(kbounds(1)-kbounds(1)+1.0d0)/log(10.0d0), log(kbounds(2)-kbounds(1)+1.0d0)/log(10.0d0), nk)
    knotsk = knotsk + (kbounds(1) - 1.0d0)
    knotsb = linspace(kbounds(1), kbounds(2), nb)

    invTk = spbas(rk,knotsk)
    call dgetrf(rk,rk,invTk,rk,IPIVk,INFO)
    call dgetri(rk,invTk,rk,IPIVk,WORKk,rk,INFO)

    print *, 'Calculating the steady state'
    call calcss(knotsk,knotsb,invTk,Ge,Pe,mue,vmatss,gmatss,mu0,evec,mpmat,ymat,nmat)
    muss = mu0
    call calcdiststat(knotsb,mu0,mpmat,ikvec)
    if (adjbias .eqv. .false.) evec = 0.0d0

    ! open(100, file="mu0.txt")
    ! open(101, file="mpmat.txt")
    ! open(102, file="ymat.txt")
    ! open(103, file="nmat.txt")
    !
    ! do ie = 1,ne
    !     do ib = 1,nb
    !
    !         write(100, *) mu0(ib,ie)
    !         write(101, *) mpmat(ib,ie)
    !         write(102, *) ymat(ib,ie)
    !         write(103, *) nmat(ib,ie)
    !
    !     end do
    ! end do
    !
    ! open(110, file="knotsb.txt")
    ! do ib = 1,nb
    !
    !     write(110, *) knotsb(ib)
    !     ! read(110, *)  knotsb(ib)
    !
    ! end do
    !
    ! close(100)
    ! close(101)
    ! close(102)
    ! close(103)
    ! close(110)

    ! fraction of e-indexed capital
    mnow = 0.0d0
    do ie = 1,ne
        psie(ie) = sum(knotsb*mu0(:,ie),1)/sum(mu0(:,ie),1)
        mnow = mnow + sum(knotsb*mu0(:,ie),1)
    end do
    psie = psie/mnow
    print *, evec, psie, mnow
    ! pause

    mbounds = (/0.75d0*mnow, 1.25d0*mnow/)
    ! mbounds = (/0.85d0, 1.65d0/)
    ! mbounds = (/1.25d0, 2.0d0/)
    knotsm = linspace(mbounds(1), mbounds(2), nm)
    invTm = spbas(rm,knotsm)
    call dgetrf(rm,rm,invTm,rm,IPIVm,INFO)
    call dgetri(rm,invTm,rm,IPIVm,WORKm,rm,INFO)

    if (fcstini) then

        call json%load_file(filename = trim(jsonfilename))
        call json%get('output.mpmat0', temp)
        mpmat0 = reshape(temp,(/nm,nz/))
        call json%get('output.pmat0', temp)
        pmat0 = reshape(temp,(/nm,nz/))

    else

        print *, 'Solving the planner''s problem'
        call solvepl(knotsk,invTk,Gz,Pz,gmatpl,nmatpl,vmatpl)
        call forecastpl(nmatpl,knotsk,knotsm,Gz,mpmat0,pmat0)

    end if

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
        !read(500, *) izvec(tt)

    end do

    ! open(500, file="izvec.txt")
    !
    ! do tt = 1,simT+drop
    !     write(500, *) izvec(tt)
    !     ! read(500, *) izvec(tt)
    ! end do
    !
    ! close(500)


    ! vmat0 = 0.0d0
    diff = 1d+4
    iter = 0

    print *, 'Solving the lumpy model'
    do while (diff>critout)

        call inner(mpmat0,pmat0,knotsk,knotsm,invTk,invTm,Gz,Pz,Ge,Pe,vmat0,gmat0,eptimein,error)

        call simulation(vmat0,mpmat0,pmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Pz,Ge,Pe,izvec,aggsim,disaggsim,mu0,eptimeout)
        Kvec = aggsim(:,6)
        Cvec = aggsim(:,4)
        call calcforecast(izvec,Kvec,1.0d0/Cvec,BetaK,Betap,R2)
        call Beta2mat(BetaK,Betap,knotsm,mpmat1,pmat1)

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
        call core%print(p, jsonfilename)
        call core%destroy(p)

    end if

    ! open(500, file="izvec.txt")
    ! open(511, file="Yvec.txt")
    ! open(512, file="Zvec.txt")
    ! open(513, file="Nvec.txt")
    ! open(514, file="Cvec.txt")
    ! open(515, file="Ivec.txt")
    ! open(516, file="Kvec.txt")
    !
    ! do tt = 1,simT+drop
    !     write(500, *) izvec(tt)
    !     write(511, *) aggsim(tt,1)
    !     write(512, *) aggsim(tt,2)
    !     write(513, *) aggsim(tt,3)
    !     write(514, *) aggsim(tt,4)
    !     write(515, *) aggsim(tt,5)
    !     write(516, *) aggsim(tt,6)
    ! end do
    !
    ! open(520, file="knotsm.txt")
    !
    ! do im = 1,nm
    !    write(520, *) knotsm(im)
    ! end do
    !
    ! open(530, file="mpmat0.txt")
    ! open(531, file="pmat0.txt")
    !
    ! do im = 1,nm
    !     do iz = 1,nz
    !         write(530, *) mpmat0(im,iz)
    !         write(531, *) pmat0(im,iz)
    !     end do
    ! end do
    !
    ! close(510)
    ! close(511)
    ! close(512)
    ! close(513)
    ! close(514)
    ! close(515)
    ! close(516)
    ! close(520)
    ! close(530)
    ! close(531)

    ! print *, 'Den Haan statistics'
    ! ! Kvec = aggsim(:,6)
    ! ! Cvec = aggsim(:,4)
    ! do tt = 1,simT+drop
    !
    !     iz = izvec(tt)
    !     ! linear interpolation
    !     ! im = gridlookup2(mnow,knotsm)
    !     ! wm = (knotsm(im+1)-mnow)/(knotsm(im+1)-knotsm(im))
    !     ! mp = wm*mpmat0(im,iz) + (1.0d0-wm)*mpmat0(im+1,iz)
    !     ! cnow = 1.0d0/(wm*pmat0(im,iz) + (1.0d0-wm)*pmat0(im+1,iz))
    !     ! log-linear interpolation
    !     im = gridlookup2(log(mnow),log(knotsm))
    !     wm = log(knotsm(im+1)/mnow)/log(knotsm(im+1)/knotsm(im))
    !     mp = exp(wm*log(mpmat0(im,iz)) + (1.0d0-wm)*log(mpmat0(im+1,iz)))
    !     cnow = 1.0d0/exp(wm*log(pmat0(im,iz)) + (1.0d0-wm)*log(pmat0(im+1,iz)))
    !
    !     mvec(tt) = mnow
    !     pvec(tt) = 1.0d0/cnow
    !     mnow = mp
    !
    ! end do
    !
    ! do iz = 1,nz
    !
    !     DHmax(iz,1)  = 100.0d0*maxval(abs(log(pack(mvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz) / pack(Kvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz) )),1)
    !     DHmean(iz,1) = 100.0d0*sum(abs(log(pack(mvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz)    / pack(Kvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz) )),1)/simT
    !     DHmax(iz,2)  = 100.0d0*maxval(abs(log(pack(pvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz) * pack(Cvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz) )),1)
    !     DHmean(iz,2) = 100.0d0*sum(abs(log(pack(pvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz)    * pack(Cvec(drop+1:simT+drop), izvec(drop+1:simT+drop)==iz) )),1)/simT
    !
    ! end do
    ! write(*,"('  max for m''  ( ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ') ')") DHmax(:,1)
    ! write(*,"('  max for p   ( ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ') ')") DHmax(:,2)
    ! write(*,"('  mean for m'' ( ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ') ')") DHmean(:,1)
    ! write(*,"('  mean for p  ( ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ', ', F8.5, ') ')") DHmean(:,2)


end program solveKS

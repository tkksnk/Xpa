program calcirf


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
    real(8) cumPz(nz,nz), rand, Zvec(irfTT), aggirf(irfTT,8), disaggirf(irfTT,7)
    ! for json
    type(json_file) :: json
    type(json_core) :: core
    real(8), allocatable :: temp(:)
    integer found
    type(json_value), pointer :: p, irf
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
    call calcss(knotsk,knotsb,invTk,Ge,Pe,mue,vmatss,gmatss,muss,evec,mpmat,ymat,nmat)
    mu0 = muss
    call calcdiststat(knotsb,mu0,mpmat,ikvec)
    if (adjbias .eqv. .false.) evec = 0.0d0

    ! fraction of e-indexed capital
    mnow = 0.0d0
    mp = 0.0d0
    ynow = 0.0d0
    nnow = 0.0d0
    do ie = 1,ne
        psie(ie) = sum(knotsb*mu0(:,ie),1)/sum(mu0(:,ie),1)
        mnow = mnow + sum(knotsb*mu0(:,ie),1)
        mp = mp + sum(mpmat(:,ie)*mu0(:,ie),1)
        ynow = ynow + sum(ymat(:,ie)*mu0(:,ie),1)
        nnow = nnow + sum(nmat(:,ie)*mu0(:,ie),1)
    end do
    psie = psie/mnow
    print *, evec, psie
    inow = GAMY*mp - (1.0d0-DELTA)*mnow
    print *, inow/mnow, mnow/ynow, nnow, 1.0d0/(ynow-inow)
    pause

    mbounds = (/0.75d0*mnow, 1.25d0*mnow/)
    knotsm = linspace(mbounds(1), mbounds(2), nm)
    invTm = spbas(rm,knotsm)
    call dgetrf(rm,rm,invTm,rm,IPIVm,INFO)
    call dgetri(rm,invTm,rm,IPIVm,WORKm,rm,INFO)

    ! print *, 'Solving the planner''s problem'
    ! call solvepl(knotsk,invTk,Gz,Pz,gmatpl,nmatpl,vmatpl)
    ! call forecastpl(nmatpl,knotsk,knotsm,Gz,mpmat0,pmat0)
    ! open(530, file="mpmat0.txt")
    ! open(531, file="pmat0.txt")
    ! read the file
    call json%initialize()
    ! call json%load_file(filename = './results_extend_KS.json')
    call json%load_file(filename = trim(jsonfilename))

    ! ! print the file to the console
    ! call json%print_file()
    ! pause

    ! extract data from the file
    ! [found can be used to check if the data was really there]
    call json%get('output.mpmat0', temp)
    mpmat0 = reshape(temp,(/nm,nz/))
    call json%get('output.pmat0', temp)
    pmat0 = reshape(temp,(/nm,nz/))

    ! open(530, file="/home/takeki/Dropbox/current/hetero4/KT180803/results/extended/KS/mpmat0.txt")
    ! open(531, file="/home/takeki/Dropbox/current/hetero4/KT180803/results/extended/KS/pmat0.txt")
    !
    ! do iz = 1,nz
    !     do im = 1,nm
    !         read(530, *) mpmat0(im,iz)
    !         read(531, *) pmat0(im,iz)
    !     end do
    ! end do
    !
    ! close(530)
    ! close(531)

    call inner(mpmat0,pmat0,knotsk,knotsm,invTk,invTm,Gz,Pz,Ge,Pe,vmat0,gmat0,eptimein,error)

    ! determine the path of TFP shock
    Zvec(1:irfdrop+1) = 1.0d0
    Zvec(irfdrop+2) = max(exp(RHO*log(Zvec(irfdrop+1))+shocksize),Gz(1))

	do tt=irfdrop+2,irfTT-1
        Zvec(tt+1) = exp(RHO*log(Zvec(tt)))
    end do

    call calcirf1(vmat0,mpmat0,pmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Pz,Ge,Pe,Zvec,aggirf,disaggirf,mu0,eptimeout)

    ! the i/o of data will be replaced by json
    ! output via json
    ! call json%create_object(p, '')
    ! how to overwrite?
    call json%remove('irf')
    call json%get(p)
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

    open(511, file="Yirvec.txt")
    open(512, file="Zirvec.txt")
    open(513, file="Nirvec.txt")
    open(514, file="Cirvec.txt")
    open(515, file="Iirvec.txt")
    open(516, file="Kirvec.txt")

    do tt = 1,irfTT
        write(511, *) aggirf(tt,1)
        write(512, *) aggirf(tt,2)
        write(513, *) aggirf(tt,3)
        write(514, *) aggirf(tt,4)
        write(515, *) aggirf(tt,5)
        write(516, *) aggirf(tt,6)
    end do

    close(511)
    close(512)
    close(513)
    close(514)
    close(515)
    close(516)


end program calcirf

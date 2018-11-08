module mod_outer


use mod_functions
use mod_spline
use mod_inner
use mod_calcss, only: calcdiststat
implicit none


contains


subroutine calcforecast(izvec,mvec,kappakp,rsq)


    use mod_parameters
    integer, intent(in) :: izvec(:)
    real(8), intent(in) :: mvec(:)
    real(8), intent(out) :: kappakp(:,:), rsq(:)
    real(8), allocatable :: XX(:,:), yy(:,:), ee(:,:), yyhat(:,:)
    real(8) yymean(1), AA(2,2), bb(2,1), det, invAA(2,2), kappa(2,1), temp(1)
    integer iz, num

    do iz = 1,nz

        num = count(izvec(drop+1:simTT-1)==iz)

        allocate(XX(num,2),yy(num,1),ee(num,1),yyhat(num,1))

        XX(:,1) = 1.0d0
        XX(:,2) = log(pack(mvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))
        yy(:,1) = log(pack(mvec(drop+2:simTT),  izvec(drop+1:simTT-1)==iz))

        yymean = sum(yy,1)/num

        AA = matmul(transpose(XX),XX)
        bb = matmul(transpose(XX),yy)
        det = AA(2,2)*AA(1,1) - AA(1,2)*AA(2,1)
        invAA(1,:) = (/AA(2,2), -AA(1,2)/)
        invAA(2,:) = (/-AA(2,1), AA(1,1)/)
        invAA =  invAA/det
        kappa = matmul(invAA,bb)

        kappakp(iz,1) = kappa(1,1)
        kappakp(iz,2) = kappa(2,1)

        yyhat(:,1) = kappa(1,1) + kappa(2,1)*XX(:,2)
        ee = yy - yyhat
        temp = matmul(reshape(ee(:,1),(/1,num/)),ee(:,1))
        temp = temp/matmul(reshape(yy(:,1)-yymean(1),(/1,num/)),yy(:,1)-yymean(1))
        rsq(iz) = 1.0d0-temp(1)

        deallocate(XX,yy,ee,yyhat)

    end do


end subroutine calcforecast


subroutine outer(vmat0,gmat0,mpmat0,knotsk,knotsm,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,psix,muxz,evec,mpmat1,eptime)


    use mod_parameters
    real(8), intent(in) :: vmat0(:,:,:,:), gmat0(:,:,:,:)
    real(8), intent(in) :: mpmat0(:,:), knotsk(:), knotsm(:), invTk(:,:), invTm(:,:), Gz(:), Gl(:), Pz(:,:), Gez(:,:), Gy(:), Gd(:), Pez(:,:,:), Py(:,:), Pd(:,:), &
        psix(:,:), muxz(:,:), evec(:,:)
    integer, intent(in) :: idmat(:,:)
    real(8), intent(out) :: mpmat1(:,:), eptime
    real(8) vcond(nk,nm), EVmate(nk,nm,nz,nx), mnow, mp, wm, vm(nk), znow, lnow, enow, ynow, beta, rr, wr, yterm, know, klow, khigh, kp, df, d2f, mpvec(nx)
    real(8), allocatable :: cfmat(:,:,:), cfmate(:,:,:,:,:)
    integer i, ix, jx, ie, je, iy, jy, id, jd, iz, im, jm, c1, c2, cr


    call system_clock(c1,cr)

    if (linflag) then
        allocate(cfmat(4,rk+1,1),cfmate(4,rk+1,nm,nz,nx))
    else
        allocate(cfmat(16,rk+1,rm+1),cfmate(16,rk+1,rm+1,nz,nx))
    end if

    ! fit splines
    if (spliflag) then

        do ix = 1,nx

            do iz = 1,nz

                do im = 1,nm

                    if (linflag) then
                        cfmate(:,:,im,iz,ix) = spfit(invTk,reshape(gmat0(:,im,iz,ix),(/nk/)),rk,knotsk)
                    else
                        cfmate(:,:,:,iz,ix) = spfit2(invTk,invTm,reshape(gmat0(:,:,iz,ix),(/nk,nm/)),rk,rm,knotsk,knotsm)
                    end if

                end do

            end do

        end do

    else

        !$omp parallel do private(ix,ie,iy,id,iz,vcond,jx,je,jy,jd)
        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            do iz = 1,nz

                vcond = 0.0d0
                do jx = 1,nx

                    je = idmat(jx,1)
                    jy = idmat(jx,2)
                    jd = idmat(jx,3)
                    vcond = vcond + Py(iy,jy)*Pd(id,jd)*(Pez(ie,je,2*(iz-1)+1)*vmat0(:,:,1,jx) + Pez(ie,je,2*(iz-1)+2)*vmat0(:,:,2,jx))

                end do

                EVmate(:,:,iz,ix) = vcond
                if (linflag .eqv. .false.) cfmate(:,:,:,iz,ix) = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

            end do

        end do

    end if

    !omp parallel do private(im,iz,znow,lnow,mnow,mp,rr,wr,ix,know,kp,df,d2f,vcond,jm,wm,vm,cfmat,ie,iy,id,enow,ynow,beta,yterm,klow,khigh)
    do im = 1,nm

        do iz = 1,nz

            znow = Gz(iz)
            lnow = Gl(iz)
            mnow = knotsm(im)
            mp = mpmat0(im,iz)
            rr = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
            wr = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

            ! evaluate k' = g(k,m,e,z) by using splines
            do ix = 1,nx

                ! explicit aggregation: aggregate capital indexed by individual productivity
                if (naiveflag) then
                    know = mnow
                else
                    know = psix(ix,iz)*mnow
                end if

                if (spliflag) then

                    if (linflag) then
                        call speva(reshape(cfmate(:,:,im,iz,ix),(/4,rk+1/)),know,rk,knotsk,kp,df,d2f)
                    else
                        call speva2(reshape(cfmate(:,:,:,iz,ix),(/16,rk+1,rm+1/)),know,mnow,rk,rm,knotsk,knotsm,kp,df,d2f)
                    end if

                else

                    ! linear interpolation
                    if (linflag) then
                        vcond = EVmate(:,:,iz,ix)
                        jm = gridlookup2(mp,knotsm)
                        wm = (knotsm(jm+1)-mp)/(knotsm(jm+1)-knotsm(jm))
                        vm = wm*vcond(:,jm) + (1.0d0-wm)*vcond(:,jm+1)
                        ! fit univeriate splines
                        cfmat(:,:,1) = spfit(invTk,vm,rk,knotsk)
                    else
                        cfmat = cfmate(:,:,:,iz,ix)
                    end if

                    ie = idmat(ix,1)
                    iy = idmat(ix,2)
                    id = idmat(ix,3)

                    enow = Gez(ie,iz)
                    ynow = Gy(iy)
                    beta = Gd(id)
                    yterm = wr*ynow*enow + rr*know/THETA

                    klow = knotsk(1)
                    khigh = min(yterm,knotsk(nk))
                    if (nraflag) then
                        call nra(know,klow,khigh,kp,mp,cfmat,yterm,knotsk,knotsm,beta)
                    else
                        call gss(klow,khigh,kp,mp,cfmat,yterm,knotsk,knotsm,beta)
                    end if

                end if

                ! end-of-period distribution
                mpvec(ix) = kp*THETA + evec(ix,iz) ! only THETA fraction will be active

            end do

            mpmat1(im,iz) = sum(muxz(:,iz)*mpvec,1)

        end do ! z

    end do ! m

    call system_clock(c2, cr)
    eptime = dble(c2-c1)/cr
    ! write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine outer


subroutine simulation(vmat0,gmat0,mpmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Gl,Pz,Gez,Gy,Gd,Pez,Py,Pd,idmat,izvec,aggsim,disaggsim,mu0,eptime)


    use mod_parameters
    real(8), intent(in) :: vmat0(:,:,:,:), gmat0(:,:,:,:)
    real(8), intent(in) :: mpmat0(:,:), knotsk(:), knotsm(:), knotsb(:), invTk(:,:), invTm(:,:), Gz(:), Gl(:), Pz(:,:), Gez(:,:), Gy(:), Gd(:), Pez(:,:,:), Py(:,:), Pd(:,:)
    integer, intent(in) :: idmat(:,:), izvec(:)
    real(8), intent(inout) :: mu0(:,:)
    real(8), intent(out) :: aggsim(:,:), disaggsim(:,:), eptime
    real(8) vcond(nk,nm), EVmate(nk,nm,nz,nx), Pe(ne,ne), mnow, mp, wm, vm(nk), &
        znow, lnow, enow, ynow, beta, r0, w0, yterm, know, klow, khigh, kp, df, d2f, kp1, kp2, wb(nb,nx), mvec(nx), mu1(nb,nx), mpmat(nb,nx)
    real(8), allocatable :: cfmat(:,:,:), cfmate(:,:,:,:,:)
    integer ix, ie, iy, id, iz, im, ib, tt, jx, je, jy, jd, jz, jm, kb(nb,nx), c1, c2, cr, tott
    integer ThresholdWvec(7)
    real(8) ShareWvec(8), gini


    call system_clock(c1,cr)

    if (linflag) then
        allocate(cfmat(4,rk+1,1),cfmate(4,rk+1,nm,nz,nx))
    else
        allocate(cfmat(16,rk+1,rm+1),cfmate(16,rk+1,rm+1,nz,nx))
    end if

    ! fit splines
    if (spliflag) then

        do ix = 1,nx

            do iz = 1,nz

                do im = 1,nm

                    if (linflag) then
                        cfmate(:,:,im,iz,ix) = spfit(invTk,reshape(gmat0(:,im,iz,ix),(/nk/)),rk,knotsk)
                    else
                        cfmate(:,:,:,iz,ix) = spfit2(invTk,invTm,reshape(gmat0(:,:,iz,ix),(/nk,nm/)),rk,rm,knotsk,knotsm)
                    end if

                end do

            end do

        end do

    else

        !$omp parallel do private(ix,ie,iy,id,iz,vcond,jx,je,jy,jd)
        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            do iz = 1,nz

                vcond = 0.0d0
                do jx = 1,nx

                    je = idmat(jx,1)
                    jy = idmat(jx,2)
                    jd = idmat(jx,3)
                    vcond = vcond + Py(iy,jy)*Pd(id,jd)*(Pez(ie,je,2*(iz-1)+1)*vmat0(:,:,1,jx) + Pez(ie,je,2*(iz-1)+2)*vmat0(:,:,2,jx))

                end do

                EVmate(:,:,iz,ix) = vcond
                if (linflag .eqv. .false.) cfmate(:,:,:,iz,ix) = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

            end do

        end do

    end if

    tott = size(izvec,1)
    do tt = 1,tott

        iz = izvec(tt)    ! index of z
        jz = izvec(tt+1)  ! index of z'
        ! transition matrix given z and z'
        if (iz==1 .and. jz==1) Pe = Pez(:,:,1)/Pz(1,1)
        if (iz==1 .and. jz==2) Pe = Pez(:,:,2)/Pz(1,2)
        if (iz==2 .and. jz==1) Pe = Pez(:,:,3)/Pz(2,1)
        if (iz==2 .and. jz==2) Pe = Pez(:,:,4)/Pz(2,2)

        znow = Gz(iz)
        lnow = Gl(iz)

        mnow = 0.0d0
        do ix = 1,nx

            mnow = mnow + sum(knotsb*mu0(:,ix),1)
            mvec(ix) = sum(knotsb*mu0(:,ix),1)/sum(mu0(:,ix),1)

        end do

        r0 = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
        w0 = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

        if (spliflag .eqv. .true. .and. linflag .eqv. .true.) then
            ! linear interpolation for mnow
            im = gridlookup2(mnow,knotsm)
            wm = (knotsm(im+1)-mnow)/(knotsm(im+1)-knotsm(im))
            mp = wm*mpmat0(im,iz) + (1.0d0-wm)*mpmat0(im+1,iz)
        else
            ! log-linear interpolation for mnow
            im = gridlookup2(log(mnow),log(knotsm))
            wm = log(knotsm(im+1)/mnow)/log(knotsm(im+1)/knotsm(im))
            mp = exp(wm*log(mpmat0(im,iz)) + (1.0d0-wm)*log(mpmat0(im+1,iz)))
        end if

        mu1 = 0.0d0
        ! NOTE: 180816 added
        mpmat = 0.0d0
        kb = 0
        wb = 0.0d0
        ! NOTE: wm is included?
        !$omp parallel do private(ix,ib,know,kp1,kp2,df,d2f,kp,vcond,jm,wm,vm,cfmat,ie,iy,id,enow,ynow,beta,yterm,klow,khigh)
        do ix = 1,nx

            do ib = 1,nb

                if (mu0(ib,ix)>0.0d0) then

                    know = knotsb(ib)

                    if (spliflag) then

                        if (linflag) then
                            call speva(reshape(cfmate(:,:,im,iz,ix),(/4,rk+1/)),know,rk,knotsk,kp1,df,d2f)
                            call speva(reshape(cfmate(:,:,im+1,iz,ix),(/4,rk+1/)),know,rk,knotsk,kp2,df,d2f)
                            kp = wm*kp1 + (1.0d0-wm)*kp2
                            ! kp = exp(wm*log(kp1) + (1.0d0-wm)*log(kp2))
                        else
                            call speva2(reshape(cfmate(:,:,:,iz,ix),(/16,rk+1,rm+1/)),know,mnow,rk,rm,knotsk,knotsm,kp,df,d2f)
                        end if

                    else

                    ! linear interpolation for mp
                    if (linflag) then
                        vcond = EVmate(:,:,iz,ix)
                        jm = gridlookup2(mp,knotsm)
                        wm = (knotsm(jm+1)-mp)/(knotsm(jm+1)-knotsm(jm))
                        vm = wm*vcond(:,jm) + (1.0d0-wm)*vcond(:,jm+1)
                        ! fit univeriate splines
                        cfmat(:,:,1) = spfit(invTk,vm,rk,knotsk)
                    else
                        cfmat = cfmate(:,:,:,iz,ix)
                    end if

                        ie = idmat(ix,1)
                        iy = idmat(ix,2)
                        id = idmat(ix,3)

                        enow = Gez(ie,iz)
                        ynow = Gy(iy)
                        beta = Gd(id)
                        yterm = w0*ynow*enow + r0*know/THETA

                        klow = knotsk(1)
                        khigh = min(yterm,knotsk(nk))
                        if (nraflag) then
                            call nra(know,klow,khigh,kp,mp,cfmat,yterm,knotsk,knotsm,beta)
                        else
                            call gss(klow,khigh,kp,mp,cfmat,yterm,knotsk,knotsm,beta)
                        end if

                    end if

                    mpmat(ib,ix) = kp*THETA ! NOTE: 180816 added
                    kb(ib,ix) = gridlookup2(kp,knotsb)
                    wb(ib,ix) = (knotsb(kb(ib,ix)+1)-kp)/(knotsb(kb(ib,ix)+1)-knotsb(kb(ib,ix)))

                end if

            end do

        end do

        ! update the disribution using k1 and w1
        mu1 = 0.0d0
        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            do ib = 1,nb

                if (mu0(ib,ix)>0.0d0) then

                    do jx = 1,nx

                        je = idmat(jx,1)
                        jy = idmat(jx,2)
                        jd = idmat(jx,3)

                        mu1(kb(ib,ix),jx)   = mu1(kb(ib,ix),jx)   + mu0(ib,ix)*wb(ib,ix)*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA
                        mu1(kb(ib,ix)+1,jx) = mu1(kb(ib,ix)+1,jx) + mu0(ib,ix)*(1.0d0-wb(ib,ix))*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA
                        ! newborns
                        mu1(1,jx) = mu1(1,jx) + mu0(ib,ix)*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*(1.0d0-THETA)

                    end do

                end if

            end do

        end do

        ! calculate aggregate variables
        aggsim(tt,2) = znow
        aggsim(tt,6) = mnow
        aggsim(tt,1) = znow*mnow**ALPHA*lnow**(1.0d0-ALPHA) !sum(sum(mu0*ymat,1),1)
        aggsim(tt,3) = lnow !sum(sum(mu0*nmat,1),1)
        aggsim(tt,4) = znow*mnow**ALPHA*lnow**(1.0d0-ALPHA)+(1.0d0-DELTA)*mnow-sum(sum(mu0*(mpmat),1),1)
        aggsim(tt,5) = sum(sum(mu0*mpmat,1),1)-(1.0d0-DELTA)*mnow !sum(sum(mu0*imat,1),1)
        aggsim(tt,7) = sum(sum(mu0*mpmat,1),1)
        aggsim(tt,8) = sum(sum(mu0*mpmat,1),1)/mnow-(1.0d0-DELTA)

        ! calculate disaggregate variables
        call calcdiststat(knotsb,mu0,ThresholdWvec,ShareWvec,gini)
        disaggsim(tt,1:8) = ShareWvec
        disaggsim(tt,9) = gini

        if (mod(tt,diagnumout)==0) then

            write(*,"('  time = ', I5, ': m =', F10.5, ' sum(mu) = ', F10.5)") tt, mnow, sum(sum(mu1,1),1)

        end if

        ! update distribution
        ! NOTE: normalization needed?
        if (spliflag) then
            mu0 = mu1/sum(sum(mu1,1),1)
        else
            mu0 = mu1
        end if

    end do ! end mu loop


    call system_clock(c2, cr)
    eptime = dble(c2-c1)/cr
    ! write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine simulation


end module mod_outer

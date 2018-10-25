module mod_simulation


use mod_functions
use mod_spline
use mod_calcss, only: calcdiststat
use mod_inner
implicit none


contains


subroutine Beta2mat(BetaK,Betap,knotsm,mpmat,pmat)


    use mod_parameters
    real(8), intent(in) :: BetaK(nz,2), Betap(nz,2), knotsm(nm)
    real(8), intent(out) :: mpmat(nm,nz), pmat(nm,nz)
    real(8) mnow
    integer iz, im


    do iz = 1,nz

        do im = 1,nm

            mnow = knotsm(im)
            mpmat(im,iz) = exp( BetaK(iz,1) + BetaK(iz,2)*log(mnow) )
            pmat(im,iz)  = exp( Betap(iz,1) + Betap(iz,2)*log(mnow) )

        end do

    end do


end subroutine Beta2mat


subroutine calcforecast(izvec,mvec,pvec,kappakp,kappap,rsq)


    use mod_parameters
    integer, intent(in) :: izvec(:)
    real(8), intent(in) :: mvec(:), pvec(:)
    real(8), intent(out) :: kappakp(:,:), kappap(:,:), rsq(:,:)
    real(8), allocatable :: XX(:,:), yy(:,:), ee(:,:), yyhat(:,:)
    real(8) yymean(2), AA(2,2), bb(2,2), det, invAA(2,2), kappa(2,2), temp(1)
    integer iz, num

    do iz = 1,nz

        num = count(izvec(drop+1:simTT-1)==iz)

        allocate(XX(num,2),yy(num,2),ee(num,2),yyhat(num,2))

        XX(:,1) = 1.0d0
        XX(:,2) = log(pack(mvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))
        yy(:,1) = log(pack(mvec(drop+2:simTT),  izvec(drop+1:simTT-1)==iz))
        yy(:,2) = log(pack(pvec(drop+1:simTT-1),izvec(drop+1:simTT-1)==iz))

        AA = matmul(transpose(XX),XX)
        bb = matmul(transpose(XX),yy)
        det = AA(2,2)*AA(1,1) - AA(1,2)*AA(2,1)
        invAA(1,:) = (/AA(2,2), -AA(1,2)/)
        invAA(2,:) = (/-AA(2,1), AA(1,1)/)
        invAA =  invAA/det
        kappa = matmul(invAA,bb)

        kappakp(iz,1) = kappa(1,1)
        kappakp(iz,2) = kappa(2,1)
        kappap(iz,1)  = kappa(1,2)
        kappap(iz,2)  = kappa(2,2)

        yymean = sum(yy,1)/num
        yyhat(:,1) = kappa(1,1) + kappa(2,1)*XX(:,2)
        yyhat(:,2) = kappa(1,2) + kappa(2,2)*XX(:,2)
        ee = yy - yyhat
        temp = matmul(reshape(ee(:,1),(/1,num/)),ee(:,1))
        temp = temp/matmul(reshape(yy(:,1)-yymean(1),(/1,num/)),yy(:,1)-yymean(1))
        rsq(iz,1) = 1.0d0-temp(1)
        temp = matmul(reshape(ee(:,2),(/1,num/)),ee(:,2))
        temp = temp/matmul(reshape(yy(:,2)-yymean(2),(/1,num/)),yy(:,2)-yymean(2))
        rsq(iz,2) = 1.0d0-temp(1)

        deallocate(XX,yy,ee,yyhat)

    end do


end subroutine calcforecast


subroutine calcirf1(vmat0,mpmat0,pmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Pz,Ge,Pe,Zvec,aggirf,disaggirf,mu0,eptime)


    use mod_parameters
    real(8), intent(in) :: vmat0(:,:,:,:), mpmat0(:,:), pmat0(:,:), knotsk(:), knotsm(:), knotsb(:), invTk(:,:), invTm(:,:), &
    Gz(:), Pz(:,:), Ge(:), Pe(:,:), Zvec(:)
    real(8), intent(out) :: aggirf(:,:), disaggirf(:,:), eptime
    real(8), intent(inout) :: mu0(:,:)
    real(8) vcond(nk,nm), cfmatz(16,rk+1,rm+1,nz,ne), cfmate(16,rk+1,rm+1,ne), cfmat(16,rk+1,rm+1)
    real(8) mnow, znow, mp, klow, khigh, ev, edv, ed2v, pL, pH, p0, p1, B0, diff, wm, plow, phigh
    integer im, iz, jz, kz, ie, je, ib, tt, iterbp, ct1, ct2, cr
    ! NOTE: the name of mpmat is confusing
    real(8) bmat(nb,ne), mpmat(nb,ne), ymat(nb,ne), imat(nb,ne), nmat(nb,ne), mu1(nb,ne)
    real(8) wz, aggirfz(2,8), disaggirfz(2,7), Yvecz(2), Nvecz(2), Cvecz(2), Ivecz(2), mu1z(nb,ne,2)

#ifdef MATLAB_MEX_FILE
    integer, external :: mexPrintf
    integer k
    character(len=80) line
    LOGICAL :: INTERRUPTED
#endif


    call system_clock(ct1, cr)

    ! grid on distribution
    do ie = 1,ne
        bmat(:,ie) = knotsb
    end do

    ! fit splines to conditional expectation
    do iz = 1,nz

        do ie = 1,ne

            vcond = 0.0d0

            do jz = 1,nz

                do je = 1,ne

                    vcond = vcond + Pz(iz,jz)*Pe(ie,je)*reshape(vmat0(:,:,jz,je),(/nk,nm/))

                end do

            end do

            cfmatz(:,:,:,iz,ie) = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

        end do

    end do

    ! simulation
    do tt = 1,irfTT

        kz = gridlookup2(Zvec(tt),Gz)
        wz = (Gz(kz+1)-Zvec(tt))/(Gz(kz+1)-Gz(kz))

		! calculate mu1(iz) & mu1(iz+1), then linearly interpolate them
        ! NOTE: maybe it is better to interpolate inside pricemap?
		do iz = kz,kz+1

            znow = Gz(iz)
            mnow = sum(sum(mu0*bmat,1),1)
            ! linear interpolation
            ! im = gridlookup2(mnow,knotsm)
            ! wm = (knotsm(im+1)-mnow)/(knotsm(im+1)-knotsm(im))
            ! mp = wm*mpmat0(im,iz) + (1.0d0-wm)*mpmat0(im+1,iz)
            ! log-linear interpolation
            im = gridlookup2(log(mnow),log(knotsm))
            wm = log(knotsm(im+1)/mnow)/log(knotsm(im+1)/knotsm(im))
            mp = exp(wm*log(mpmat0(im,iz)) + (1.0d0-wm)*log(mpmat0(im+1,iz)))
            cfmate = reshape(cfmatz(:,:,:,iz,:), (/16,rk+1,rm+1,ne/))

            if (outermkt .eqv. .false.) then

                p0 = exp(wm*log(pmat0(im,iz)) + (1.0d0-wm)*log(pmat0(im+1,iz)))
                call pricemapks(p0,knotsk,knotsm,Ge,Pe,znow,mnow,mp,knotsb,mu0,cfmate,p1,mpmat,ymat,imat,nmat,mu1)

            else

                ! bisection over p NOTE: if instead w?
                cfmat = reshape(cfmate(:,:,:,ceiling(dble(ne)/2.0d0)), (/16,rk+1,rm+1/))
                klow = 0.5d0*mp ! low k -> high edv
                call speva2(cfmat,klow,mp,rk,rm,knotsk,knotsm,ev,edv,ed2v)
                pH = BETA*edv/GAMY ! FOC of k

                khigh = 1.5d0*mp ! high k -> low edv
                call speva2(cfmat,khigh,mp,rk,rm,knotsk,knotsm,ev,edv,ed2v)
                pL = BETA*edv/GAMY

                ! preserve the initial values for display
                plow = pL
                phigh = pH

                diff = 1d+4
                iterbp = 0

                do while (diff>critbp .and. iterbp<50)

                    p0 = (pL+pH)/2.0d0
                    call pricemapks(p0,knotsk,knotsm,Ge,Pe,znow,mnow,mp,knotsb,mu0,cfmate,p1,mpmat,ymat,imat,nmat,mu1)

                    B0 = p0-p1

                    if (B0<0) then
                        pL = p0
                    else
                        pH = p0
                    end if

                    diff = pH-pL
                    iterbp = iterbp + 1
                    ! diagnosis
                    ! write(*,"('  bisection ', I4, '  pH-pL = ', F10.5)") iterbp, diff
                    ! print *, p0

                end do

            end if

            ! calculate aggregate variables
            aggirfz(iz-kz+1,2) = znow
            aggirfz(iz-kz+1,6) = mnow
            aggirfz(iz-kz+1,1) = sum(sum(mu0*ymat,1),1)
            aggirfz(iz-kz+1,3) = sum(sum(mu0*nmat,1),1)
            aggirfz(iz-kz+1,4) = sum(sum(mu0*(ymat-imat),1),1)
            aggirfz(iz-kz+1,5) = sum(sum(mu0*imat,1),1)
            aggirfz(iz-kz+1,7) = sum(sum(mu0*mpmat,1),1)
            aggirfz(iz-kz+1,8) = sum(sum(mu0*imat,1),1)/mnow
            ! ! calculate aggregate variables
            ! Yvecz(iz-kz+1) = sum(sum(mu0*ymat,1),1)
            ! Nvecz(iz-kz+1) = sum(sum(mu0*nmat,1),1)
            ! Cvecz(iz-kz+1) = sum(sum(mu0*(ymat-imat),1),1)
            ! Ivecz(iz-kz+1) = sum(sum(mu0*imat,1),1)
            mu1z(:,:,iz-kz+1) = mu1

            ! calculate disaggregate variables
            call calcdiststat(knotsb,mu0,mpmat,disaggirfz(iz-kz+1,:))

            ! diagnosis
            ! if (mod(tt,1)==0) then
            if (mod(tt,1)==0) then
                write(*,"('  at iz = ', I5, ': pnow =', F8.5, ', pl =', F8.5, ', ph =', F8.5, ' ', I2)") iz, p1, plow, phigh, iterbp
                ! write(*,"(F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)") irfmat(tt,:)
            end if

        end do ! iz

        ! linear interpolation
        aggirf(tt,:) = wz*aggirfz(1,:) + (1.0d0-wz)*aggirfz(2,:)
        disaggirf(tt,:) = wz*disaggirfz(1,:) + (1.0d0-wz)*disaggirfz(2,:)
        ! aggirf(tt,2) = wz*Gz(kz) + (1.0d0-wz)*Gz(kz+1)
        ! aggirf(tt,6) = mnow
        ! aggirf(tt,1) = wz*Yvecz(1) + (1.0d0-wz)*Yvecz(2)
        ! aggirf(tt,3) = wz*Nvecz(1) + (1.0d0-wz)*Nvecz(2)
        ! aggirf(tt,4) = wz*Cvecz(1) + (1.0d0-wz)*Cvecz(2)
        ! aggirf(tt,5) = wz*Ivecz(1) + (1.0d0-wz)*Ivecz(2)

        ! diagnosis
        ! if (mod(tt,1)==0) then
        if (mod(tt,1)==0) then
            ! write(*,"('  time = ', I5, ': pnow =', F8.5, ', pl =', F8.5, ', ph =', F8.5, ' ', I2)") tt, p1, plow, phigh, iterbp
            write(*,"('  time = ', I5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)") tt, aggirf(tt,1:6)
            print *, disaggirf(tt,:)
        end if

        ! update distribution
        mu0 = wz*mu1z(:,:,1) + (1.0d0-wz)*mu1z(:,:,2)

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine calcirf1


subroutine simulation(vmat0,mpmat0,pmat0,knotsk,knotsm,knotsb,invTk,invTm,Gz,Pz,Ge,Pe,izvec,aggsim,disaggsim,mu0,eptime)


    use mod_parameters
    real(8), intent(in) :: vmat0(:,:,:,:), mpmat0(:,:), pmat0(:,:), knotsk(:), knotsm(:), knotsb(:), invTk(:,:), invTm(:,:), &
    Gz(:), Pz(:,:), Ge(:), Pe(:,:)
    integer, intent(in) :: izvec(:)
    real(8), intent(out) :: aggsim(:,:), disaggsim(:,:), eptime
    real(8), intent(inout) :: mu0(:,:)
    real(8) vcond(nk,nm), cfmatz(16,rk+1,rm+1,nz,ne), cfmate(16,rk+1,rm+1,ne), cfmat(16,rk+1,rm+1)
    real(8) mnow, znow, mp, klow, khigh, ev, edv, ed2v, pL, pH, p0, p1, B0, diff, wm, plow, phigh
    integer im, iz, jz, ie, je, ib, tt, iterbp, ct1, ct2, cr
    ! NOTE: the name of mpmat is confusing
    real(8) bmat(nb,ne), mpmat(nb,ne), ymat(nb,ne), imat(nb,ne), nmat(nb,ne), mu1(nb,ne)

#ifdef MATLAB_MEX_FILE
    integer, external :: mexPrintf
    integer k
    character(len=80) line
    LOGICAL :: INTERRUPTED
#endif


    call system_clock(ct1, cr)

    ! grid on distribution
    do ie = 1,ne
        bmat(:,ie) = knotsb
    end do

    ! fit splines to conditional expectation
    do iz = 1,nz

        do ie = 1,ne

            vcond = 0.0d0

            do jz = 1,nz

                do je = 1,ne

                    vcond = vcond + Pz(iz,jz)*Pe(ie,je)*reshape(vmat0(:,:,jz,je),(/nk,nm/))

                end do

            end do

            cfmatz(:,:,:,iz,ie) = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

        end do

    end do

    ! simulation
    do tt = 1,simTT

        iz = izvec(tt)
        znow = Gz(iz)
        mnow = sum(sum(mu0*bmat,1),1)
        ! linear interpolation
        ! im = gridlookup2(mnow,knotsm)
        ! wm = (knotsm(im+1)-mnow)/(knotsm(im+1)-knotsm(im))
        ! mp = wm*mpmat0(im,iz) + (1.0d0-wm)*mpmat0(im+1,iz)
        ! log-linear interpolation
        im = gridlookup2(log(mnow),log(knotsm))
        wm = log(knotsm(im+1)/mnow)/log(knotsm(im+1)/knotsm(im))
        mp = exp(wm*log(mpmat0(im,iz)) + (1.0d0-wm)*log(mpmat0(im+1,iz)))
        cfmate = reshape(cfmatz(:,:,:,iz,:), (/16,rk+1,rm+1,ne/))

        if (outermkt .eqv. .false.) then

            p0 = exp(wm*log(pmat0(im,iz)) + (1.0d0-wm)*log(pmat0(im+1,iz)))
            call pricemapks(p0,knotsk,knotsm,Ge,Pe,znow,mnow,mp,knotsb,mu0,cfmate,p1,mpmat,ymat,imat,nmat,mu1)

        else

            ! bisection over p NOTE: if instead w?
            cfmat = reshape(cfmate(:,:,:,ceiling(dble(ne)/2.0d0)), (/16,rk+1,rm+1/))
            klow = 0.5d0*mp ! low k -> high edv
            call speva2(cfmat,klow,mp,rk,rm,knotsk,knotsm,ev,edv,ed2v)
            pH = BETA*edv/GAMY ! FOC of k

            khigh = 1.5d0*mp ! high k -> low edv
            call speva2(cfmat,khigh,mp,rk,rm,knotsk,knotsm,ev,edv,ed2v)
            pL = BETA*edv/GAMY

            ! preserve the initial values for display
            plow = pL
            phigh = pH

            diff = 1d+4
            iterbp = 0

            do while (diff>critbp .and. iterbp<50)

                p0 = (pL+pH)/2.0d0
                call pricemapks(p0,knotsk,knotsm,Ge,Pe,znow,mnow,mp,knotsb,mu0,cfmate,p1,mpmat,ymat,imat,nmat,mu1)

                B0 = p0-p1

                if (B0<0) then
                    pL = p0
                else
                    pH = p0
                end if

                diff = pH-pL
                iterbp = iterbp + 1
                ! diagnosis
                ! write(*,"('  bisection ', I4, '  pH-pL = ', F10.5)") iterbp, diff
                ! print *, p0

            end do

        end if

        ! calculate aggregate variables
        aggsim(tt,2) = znow
        aggsim(tt,6) = mnow
        aggsim(tt,1) = sum(sum(mu0*ymat,1),1)
        aggsim(tt,3) = sum(sum(mu0*nmat,1),1)
        aggsim(tt,4) = sum(sum(mu0*(ymat-imat),1),1)
        aggsim(tt,5) = sum(sum(mu0*imat,1),1)
        aggsim(tt,7) = sum(sum(mu0*mpmat,1),1)
        aggsim(tt,8) = sum(sum(mu0*imat,1),1)/mnow

        ! calculate disaggregate variables
        call calcdiststat(knotsb,mu0,mpmat,disaggsim(tt,:))

        ! diagnosis
        ! if (mod(tt,1)==0) then
        if (mod(tt,diagnumout)==0) then
#ifdef MATLAB_MEX_FILE
            write(line,"('  time      ', I4, '  pnow =', F8.5, ', pl =', F8.5, ', ph =', F8.5, ' ', I2)") tt, p1, plow, phigh, iterbp
            k = mexPrintf(line//achar(10))
            ! write(line,"(F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)") aggmat(tt,:)
    		k = mexPrintf(line//achar(10))
    		call mexEvalString("drawnow")
#else
            write(*,"('  time = ', I5, ': pnow =', F8.5, ', pl =', F8.5, ', ph =', F8.5, ' ', I2)") tt, p1, plow, phigh, iterbp
            ! print *, mumat(tt,:)
            print *, sum(sum(mu0,1),1)
            ! write(*,"(F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)") aggmat(tt,:)
#endif
        end if

#ifdef MATLAB_MEX_FILE
		INTERRUPTED = utIsInterruptPendingInFortran()
		IF (INTERRUPTED) RETURN
#endif

        ! update distribution
        mu0 = mu1

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine simulation


subroutine pricemapks(p0,knotsk,knotsm,Ge,Pe,znow,mnow,mp,knotsb,mu0,cfmate,p1,mpmat,ymat,imat,nmat,mu1)


    use mod_parameters
    real(8), intent(in) :: p0, knotsk(:), knotsm(:), Ge(ne), Pe(ne,ne), znow, mnow, mp, knotsb(:), mu0(:,:), cfmate(:,:,:,:)
    real(8), intent(out) :: p1, mpmat(:,:), ymat(:,:), imat(:,:), nmat(:,:), mu1(:,:)
!    integer, intent(out) :: error
    real(8) cfmat(16,rk+1,rm+1)
    real(8) w0, enow, know, yterm, nnow, ynow, inow, klow, khigh
    real(8) kwnow, f0, df0, d2f0, e0now, kcnow, f1, df1, d2f1, e1now, xinow, a1(nb,ne), wb1(ne), wb2(nb,ne)
    integer kb1(ne), kb2(nb,ne)
    integer ie, ib, je


    w0 = ETA/p0
!    error = 0
    mpmat = 0.0d0 ! NOTE: added on 18/08/12
    ymat = 0.0d0
    imat = 0.0d0
    nmat = 0.0d0
    mu1 = 0.0d0

    !$omp parallel do private(enow,cfmat,kwnow,f0,df0,d2f0,e0now,know,klow,khigh,kcnow,f1,df1,d2f1,e1now,xinow,yterm,ynow,inow)
    do ie = 1,ne

        enow = Ge(ie)
        cfmat = reshape(cfmate(:,:,:,ie), (/16,rk+1,rm+1/))

        if (nraflag) then
            call nra(mnow,knotsk(1),knotsk(nk),kwnow,cfmat,knotsk,knotsm,mp,p0)
        else
            call gss(knotsk(1),knotsk(nk),kwnow,cfmat,knotsk,knotsm,mp,p0)
        end if

        call vfuncsp2(kwnow,cfmat,knotsk,knotsm,mp,p0,f0,df0,d2f0)
        e0now = -f0 !-vfuncsp2(kwnow,cfmat,knotsk,knotsm,mp,p0)

        kb1(ie) = gridlookup2(kwnow,knotsb)
        wb1(ie) = (knotsb(kb1(ie)+1)-kwnow)/(knotsb(kb1(ie)+1)-knotsb(kb1(ie)))

        do ib = 1,nb

            know = knotsb(ib)

            if (mu0(ib,ie)>0.0d0) then

                if (B==0.0d0) then

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

                call vfuncsp2(kcnow,cfmat,knotsk,knotsm,mp,p0,f1,df1,d2f1)
                e1now = -f1 !-vfuncsp2(kcnow,cfmat,knotsk,knotsm,mp,p0)
                xinow = min( XIBAR, max(0.0d0,(e0now-e1now)/ETA) )
                a1(ib,ie) = xinow/XIBAR

                yterm = znow*enow*know**THETA
                nnow = (NU*yterm/w0)**(1.0d0/(1.0d0-NU))
                ynow = yterm*nnow**NU
                inow = GAMY*(a1(ib,ie)*kwnow + (1.0d0-a1(ib,ie))*kcnow) - (1.0d0-DELTA)*know

                ! distribution
                mpmat(ib,ie) = a1(ib,ie)*kwnow + (1.0d0-a1(ib,ie))*kcnow
                nmat(ib,ie) = nnow + xinow**2/XIBAR/2.0d0 ! NOTE: Is the adj. cost term correct?
                ymat(ib,ie) = ynow
                imat(ib,ie) = inow
                ! ikmat(ib,ie) = inow/know

                kb2(ib,ie) = gridlookup2(kcnow,knotsb)
                wb2(ib,ie) = (knotsb(kb2(ib,ie)+1)-kcnow)/(knotsb(kb2(ib,ie)+1)-knotsb(kb2(ib,ie)))

            end if

        end do

    end do
    !$omp end parallel do

    ! update the disribution using k1, w1, k2, w2, and a1
    do ie = 1,ne

        do ib = 1,nb

            if (mu0(ib,ie)>0.0d0) then

                do je = 1,ne

                    ! if adjust with prob a1(ib,ie)
                    mu1(kb1(ie),je)      = mu1(kb1(ie),je)      + mu0(ib,ie)*wb1(ie)*a1(ib,ie)*Pe(ie,je)
                    mu1(kb1(ie)+1,je)    = mu1(kb1(ie)+1,je)    + mu0(ib,ie)*(1.0d0-wb1(ie))*a1(ib,ie)*Pe(ie,je)
                    ! if not adjust with prob 1-a1(ib,ie)
                    mu1(kb2(ib,ie),je)   = mu1(kb2(ib,ie),je)   + mu0(ib,ie)*wb2(ib,ie)*(1.0d0-a1(ib,ie))*Pe(ie,je)
                    mu1(kb2(ib,ie)+1,je) = mu1(kb2(ib,ie)+1,je) + mu0(ib,ie)*(1.0d0-wb2(ib,ie))*(1.0d0-a1(ib,ie))*Pe(ie,je)

                end do

            end if

        end do

    end do

    p1 = 1.0d0/sum(sum(mu0*(ymat-imat),1),1)


end subroutine pricemapks


end module mod_simulation

module mod_outer


use mod_spline
use mod_inner
implicit none


contains


subroutine outer(vmat0,mpmat0,knotsk,knotsm,invTk,invTm,Gz,Pz,Ge,Pe,psie,mue,evec,mpmat1,pmat1,eptime)


    use mod_parameters
    real(8), intent(in)  :: vmat0(:,:,:,:), mpmat0(:,:), knotsk(:), knotsm(:), invTk(:,:), invTm(:,:), &
    Gz(:), Pz(:,:), Ge(:), Pe(:,:), psie(:), mue(:), evec(:,:)
    real(8), intent(out) :: mpmat1(nm,nz), pmat1(nm,nz), eptime
    real(8) vcond(nk,nm), cfmatz(16,rk+1,rm+1,ne,nz), cfmate(16,rk+1,rm+1,ne), cfmat(16,rk+1,rm+1)
    real(8) mnow, znow, mp, klow, khigh, ev, edv, ed2v, pL, pH, p0, p1, B0, diff
    real(8) mpvec(ne), yvec(ne), Kagg, Kpagg, Yagg, Cagg, Iagg
    integer im, iz, jz, ie, je, iterbp, ct1, ct2, cr, error


    call system_clock(ct1, cr)

    ! fit splines to conditional expectation
    do iz = 1,nz

        do ie = 1,ne

            vcond = 0.0d0

            do jz = 1,nz

                do je = 1,ne

                    vcond = vcond + Pz(iz,jz)*Pe(ie,je)*reshape(vmat0(:,:,jz,je),(/nk,nm/))

                end do

            end do

            cfmatz(:,:,:,ie,iz) = spfit2(invTk,invTm,vcond,rk,rm,knotsk,knotsm)

        end do

    end do

    ! explicit aggregation of the induvidual decision rules
    mpmat1 = 0.0d0
    pmat1  = 0.0d0

    !omp parallel do private(znow,mnow,mp,cfmate,cfmat,klow,khigh,pL,pH,ev,edv,ed2v, &
    !omp diff,iterbp,p0,p1,mpvec,yvec,B0,Kagg,Kpagg,Yagg,Iagg,Cagg)
    do im = 1,nm

        do iz = 1,nz

            znow = Gz(iz)
            mnow = knotsm(im)
            mp = mpmat0(im,iz)
            cfmate = reshape(cfmatz(:,:,:,:,iz), (/16,rk+1,rm+1,ne/))

            ! bisection over p NOTE: if instead w?
            cfmat = reshape(cfmate(:,:,:,ceiling(dble(ne)/2.0d0)), (/16,rk+1,rm+1/))
            klow = 0.1d0*mp ! low k -> high edv
            call speva2(cfmat,klow,mp,rk,rm,knotsk,knotsm,ev,edv,ed2v)
            pH = BETA*edv/GAMY ! FOC of k

            khigh = 1.5d0*mp ! high k -> low edv
            call speva2(cfmat,khigh,mp,rk,rm,knotsk,knotsm,ev,edv,ed2v)
            pL = BETA*edv/GAMY

            diff = 1d+4
            iterbp = 1
            error = 0

            do while (diff>critbp .and. iterbp<50)

                p0 = (pL+pH)/2.0d0
                call pricemap(p0,znow,mnow,mp,knotsk,knotsm,Ge,Pe,psie,mue,evec,cfmate,p1,mpvec,yvec)

                B0 = p0-p1

                if (B0<0.0d0) then
                    pL = p0
                else
                    pH = p0
                end if

                diff = pH-pL
                iterbp = iterbp + 1
                ! diagnosis
                ! write(*,"('  bisection ', I4, '  pH-pL = ', F10.5)") iterbp, diff
                ! print *, p0, p1

            end do

            if (iterbp>50) error = 1

            ! NOTE: what is the relationshop between market clearing and bias correction?
            ! bias correction after market clearing?
            ! mpvec = mpvec + evec(:,1)
            ! yvec = yvec + evec(:,2)

            ! calculate aggregate variables
            Kagg = mnow
            Kpagg = sum(mpvec*mue,1)
            Yagg = sum(yvec*mue,1)
            Iagg = GAMY*Kpagg - (1.0d0-DELTA)*Kagg
            Cagg = Yagg-Iagg

            ! diagnosis
            ! write(*,"('  at (z,m) = (', F8.5, ',', F8.5, ') : pnow = ', F8.5, ', pl = ', F8.5)") znow, mnow, p1, plow

            ! NOTE: weighted average over e
            ! BUG: these two yield different results?
            mpmat1(im,iz) = Kpagg
            pmat1(im,iz)  = 1.0d0/Cagg
            ! mpmat1(im,iz) = Kpagg + sum(evec(:,1)*mue,1)
            ! pmat1(im,iz)  = 1.0d0/(Cagg + sum(evec(:,2)*mue,1))
            ! mpmat1(im,iz) = sum(mpvec*mue,1)
            ! pmat1(im,iz)  = 1.0d0/(sum((yvec-GAMY*mpvec)*mue,1)+(1.0d0-DELTA)*mnow)

        end do

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine outer


subroutine pricemap(p0,znow,mnow,mp,knotsk,knotsm,Ge,Pe,psie,mue,evec,cfmate,p1,mpvec,yvec)


    use mod_parameters
    real(8), intent(in) :: p0, znow, mnow, mp, knotsk(:), knotsm(:), Ge(:), Pe(:,:), psie(:), mue(:), evec(:,:), cfmate(:,:,:,:)
    real(8), intent(out) :: p1, mpvec(:), yvec(:)
    real(8) w0, enow, know, yterm, nnow, ynow, inow, klow, khigh
    real(8) kwnow, f0, df0, d2f0, e0now, kcnow, f1, df1, d2f1, e1now, xinow, alpha, Kpagg, Nagg, Yagg, Iagg, Cagg
    real(8) nvec(ne), ivec(ne), cfmat(16,rk+1,rm+1)
    integer ie


    w0 = ETA/p0

    !omp parallel do private(enow,cfmat,kwnow,e0now,know,klow,khigh,kcnow,e1now,xinow,alpha,yterm,nnow,ynow,inow)
    do ie = 1,ne

        ! explicit aggregation: aggregate capital indexed by individual productivity
        if (naiveflag) then
            know = mnow
        else
            know = psie(ie)*mnow
        end if

        enow = Ge(ie)
        cfmat = reshape(cfmate(:,:,:,ie),(/16,rk+1,rm+1/))

        if (nraflag) then
            call nra(mnow,knotsk(1),knotsk(nk),kwnow,cfmat,knotsk,knotsm,mp,p0)
        else
            call gss(knotsk(1),knotsk(nk),kwnow,cfmat,knotsk,knotsm,mp,p0)
        end if

        call vfuncsp2(kwnow,cfmat,knotsk,knotsm,mp,p0,f0,df0,d2f0)
        e0now = -f0

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
        e1now = -f1

        xinow = min(XIBAR, max(0.0d0,(e0now-e1now)/ETA))
        alpha = xinow/XIBAR

        yterm = znow*enow*know**THETA
        nnow = (NU*yterm/w0)**(1.0d0/(1.0d0-NU))
        ynow = yterm*nnow**NU
        inow = GAMY*(alpha*kwnow + (1.0d0-alpha)*kcnow) - (1.0d0-DELTA)*know

        mpvec(ie) = alpha*kwnow + (1.0d0-alpha)*kcnow
        yvec(ie) = ynow
        ivec(ie) = inow
        nvec(ie) = nnow + xinow**2/XIBAR/2.0d0 ! NOTE: Is the adj. cost term correct?

    end do
    !omp end parallel do

    ! NOTE: bias correction during market clearing
    if (naiveflag) then
        mpvec = mpvec + sum(mue*evec(:,1),1)
        yvec = yvec + sum(mue*evec(:,2),1)
    else
        mpvec = mpvec + evec(:,1)
        yvec = yvec + evec(:,2)
    end if

    Kpagg = sum((mpvec*mue),1)
    Yagg = sum((yvec*mue),1)
    Iagg = GAMY*Kpagg - (1.0d0-DELTA)*mnow
    Cagg = Yagg-Iagg
    p1 = 1.0d0/Cagg


end subroutine pricemap


end module mod_outer

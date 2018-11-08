module mod_calcss


use mod_functions
use mod_spline
use mkl_spblas
implicit none


contains


subroutine calcss(knotsk,knotsb,invTk,znow,lnow,Ge,Gy,Gd,Pe,Py,Pd,idmat,muxz,vmat0,gmat0,mu0,evec,mpmat)


    use mod_parameters
    real(8), intent(in) :: knotsk(:), knotsb(:), invTk(:,:), znow, lnow, Ge(:), Gy(:), Gd(:), Pe(:,:), Py(:,:), Pd(:,:), muxz(:,:)
    integer, intent(in) :: idmat(:,:)
    real(8), intent(out) :: vmat0(:,:), gmat0(:,:), mu0(:,:), evec(:,:), mpmat(:,:)
    real(8) Px(nx,nx)
    integer ik, ib, ix, jx, ie, je, iy, jy, id, jd, iter, iz, ct1, ct2, cr, i
    real(8) wb, mux(nx), muinit(nb,nx), m0, m1, diff, mL, mLnew, mH, mHnew, BL, BH, mnew, B0, rrnew
    real(8) vcond(nk), cf(4,rk+1), cfmate(4,rk+1,nx)
    real(8) beta, mnow, r0, w0, enow, ynow, know, yterm, klow, khigh, mp, df0, d2f0, mpvec(nx), eptime


    call system_clock(ct1,cr)

    ! initial aggregate capital and factor prices
    if (nd==1) then
        ! KS
        BETA = 0.989975d0
    else
        BETA = dbar
    end if

    mnow = lnow*(znow*ALPHA*BETA*THETA/(1.0d0-BETA*THETA*(1.0d0-DELTA)))**(1.0d0/(1.0d0-ALPHA))
    r0 = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
    w0 = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

    ! inital value function
    ! initial guess for vmat0 and gmat0
    do ix = 1,nx

        do ik = 1,nk

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            enow = Ge(ie)
            ynow = Gy(iy)
            know = knotsk(ik)
            yterm = w0*ynow*enow + r0*know/THETA
            vmat0(ik,ix) = log(yterm)
            gmat0(ik,ix) = 0.0d0

        end do

    end do

    ! initial distribution
    ! mu0 = 1.0d0/dble(nb*nx)
    ! NOTE: 18/08/13 initial distribution should not include high a's; otherwise they contaminate the results
    ! NOTE: Also, in each iteration (either bisection or iterative method) the initial distribution should be reset
    ! NOTE: 18/08/14 eivenvector method doesn't work as the transition matrix has no stable eivenvector...
    ! ! iterative method works with some loose criteria (say, 1e-5)
    ! ib = gridlookup2(mnow,knotsb)
    ! wb = (knotsb(ib+1)-mnow)/(knotsb(ib+1)-knotsb(ib))
    ! ! print *, ib, wb, knotsb(ib)
    ! ! pause
    mu0 = 0.0d0
    mux = (1.0d0-FractionZb)*muxz(:,1) + FractionZb*muxz(:,2)
    do ix = 1,nx
        ! NOTE: 22 Jan 2018 The initial distribution below is much more stable???
        mu0(:,ix)    = mux(ix)/dble(nb)
        ! mu0(ib,ix)   = wb*mux(ix)
        ! mu0(ib+1,ix) = (1.0d0-wb)*mux(ix)
    end do
    ! muinit = mu0

    if (bisectm) then

    	! bisection method
    	mL = 0.95d0*mnow
        ! mu0 = muinit
        call driverm(mL,znow,lnow,knotsk,knotsb,invTk,Ge,Gy,Gd,Pe,Py,Pd,idmat,0,mLnew,mpmat,vmat0,gmat0,mu0)
    	BL = mL-mLnew
        print *, BL

        ! mH = 1.05d0*mnow
    	mH = 1.15d0*mnow
        ! mH = 2.0d0*mnow
        ! mu0 = muinit
        call driverm(mH,znow,lnow,knotsk,knotsb,invTk,Ge,Gy,Gd,Pe,Py,Pd,idmat,0,mHnew,mpmat,vmat0,gmat0,mu0)
    	BH = mH-mHnew
        print *, BH

        m0 = (mL+mH)/2.0d0
    	diff = 1e+4
    	iter = 1

    	do while(diff>critbp)

            ! mu0 = muinit
            call driverm(m0,znow,lnow,knotsk,knotsb,invTk,Ge,Gy,Gd,Pe,Py,Pd,idmat,iter,m1,mpmat,vmat0,gmat0,mu0)
            ! print *, m0, m1
    	    B0 = m0-m1

    	    if (B0*BL>0.0d0) then ! BL>0 and B0>0
                mL = m0
    	        BL = B0
    	    else ! B0<0
                mH = m0
    	    end if

            diff = abs(log(m1)-log(m0))

            ! write(*,"('  bisection ', I3, '  mH-mL = ', F10.5)") iter, diff
            write(*,"('  bisection ', I3, '  diff = ', F10.5, '  m0 = ', F10.5, '  m1 = ', F10.5)") iter, diff, m0, m1
            iter = iter + 1
            ! updating price
            m0 = (mL+mH)/2.0d0

    	end do

    else
        ! iterative method
        m0 = mnow
    	diff = 1e+4
    	iter = 1

        do while (diff>critbp)

            call driverm(m0,znow,lnow,knotsk,knotsb,invTk,Ge,Gy,Gd,Pe,Py,Pd,idmat,iter,m1,mpmat,vmat0,gmat0,mu0)

            diff = abs(log(m1)-log(m0))

            write(*,"('  iteration ', I3, '  diff = ', F10.5, '  m0 = ', F10.5, '  m1 = ', F10.5)") iter, diff, m0, m1
            iter = iter + 1
            ! updating price
            ! if (diff>0.01d0) then
                m0 = dampss*m1 + (1.0d0-dampss)*m0
            ! else
            !     m0 = m1
            ! end if

        end do

    end if

    print *, 'Converged. Checking the result...'
    ! mu0 = muinit
    call driverm(m0,znow,lnow,knotsk,knotsb,invTk,Ge,Gy,Gd,Pe,Py,Pd,idmat,0,m1,mpmat,vmat0,gmat0,mu0)
    print *, abs(log(m1)-log(m0))
    ! pause

    ! calculate the bias correction terms
    mnow = m0
    r0 = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
    w0 = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

    do iz = 1,nz

        do ix = 1,nx

            ! explicit aggregation: aggregate capital indexed by individual productivity
            if (naiveflag) then
                know = mnow
            else
                know = sum(knotsb*mu0(:,ix),1)/sum(mu0(:,ix),1)
                ! know_xpa adjusted
                ! know = sum(knotsb*mu0(:,ix),1)/muxz(ix,iz)
            end if

            if (spliflag) then

                cf = spfit(invTk,gmat0(:,ix),rk,knotsk)
                call speva(cf,know,rk,knotsk,mp,df0,d2f0)

            else

                ie = idmat(ix,1)
                iy = idmat(ix,2)
                id = idmat(ix,3)

                ! compute the conditional expectaions
                vcond = 0.0d0
                do jx = 1,nx

                    je = idmat(jx,1)
                    jy = idmat(jx,2)
                    jd = idmat(jx,3)
                    vcond = vcond + Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*vmat0(:,jx)

                end do

                cf = spfit(invTk,vcond,rk,knotsk)

                ie = idmat(ix,1)
                iy = idmat(ix,2)
                id = idmat(ix,3)

                enow = Ge(ie)
                ynow = Gy(iy)
                beta = Gd(id)
                yterm = w0*ynow*enow + r0*know/THETA

                klow = knotsk(1)
                khigh = min(yterm,knotsk(nk))
                if (nraflag) then
                    call nra(know,klow,khigh,mp,cf,yterm,knotsk,beta)
                else
                    call gss(klow,khigh,mp,cf,yterm,knotsk,beta)
                end if

            end if

            mpvec(ix) = mp*THETA

        end do

        ! calculate the Jensen's inequality
        do ix = 1,nx

            evec(ix,iz) = -(mpvec(ix)-sum(mpmat(:,ix)*mu0(:,ix),1)/sum(mu0(:,ix),1))
            ! NOTE: kp adjusted, if we use the below, mvec becomes much less volatile?
            ! evec(ix,iz) = -(mpvec(ix)-sum(mpmat(:,ix)*mu0(:,ix),1)/muxz(ix,iz))

        end do

    end do

    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine calcss


subroutine driverm(m0,znow,lnow,knotsk,knotsb,invTk,Ge,Gy,Gd,Pe,Py,Pd,idmat,iterout,m1,mpmat,vmat0,gmat0,mu0)


    use mod_parameters
    real(8), intent(in) :: m0, znow, lnow, knotsk(:), knotsb(:), invTk(:,:), Ge(:), Gy(:), Gd(:), Pe(:,:), Py(:,:), Pd(:,:)
    integer, intent(in) :: idmat(:,:), iterout
    real(8), intent(inout) :: vmat0(:,:), gmat0(:,:), mu0(:,:)
    real(8), intent(out) :: m1, mpmat(:,:)
    real(8) vmat1(nk,nx), gmat1(nk,nx), vcond(nk), cf(4,rk+1), cfmate(4,rk+1,nx)
    real(8) mnow, w0, r0, enow, ynow, beta, know, yterm, klow, khigh, kp, f0, df0, d2f0
    real(8) crit, diff, diffg, diffv, eptime
    integer iter, s1, ik, ix, jx, ie, je, iy, jy, id, jd, ib, jb, kb(nb,nx), ct1, ct2, cr
    integer SOLVE(nb,nx)
    real(8) wb(nb,nx), G(nb,nx), mu1(nb,nx), temp(nb*nx)
    real(8), allocatable :: aa(:,:)
    ! for mkl_spblas
    integer, allocatable :: csrColInd(:)
    integer csrRowPtr(nb*nx+1)
    real(8), allocatable :: csrVal(:)
    !   Matrix descriptor
    TYPE(MATRIX_DESCR) descrA     ! Sparse matrix descriptor
    !   CSR matrix representation
    TYPE(SPARSE_MATRIX_T) csrA    ! Structure with sparse matrix
    integer index, i, info


    mnow = m0
    r0 = znow*lnow**(1.0d0-ALPHA)*ALPHA*mnow**(ALPHA-1.0d0) + 1.0d0 - DELTA
    w0 = znow*lnow**(-ALPHA)*(1.0d0-ALPHA)*mnow**ALPHA

    write(*,"('  iter = ', I3, '  rr = ', F10.5, '  wr = ', F10.5, '  mnow = ', F10.5)") iterout, r0-(1.0d0-DELTA), w0, m0

    vmat1 = 0.0d0
    gmat1 = 0.0d0

    diff = 1d+4
    iter = 0
    s1   = 0

    do while (diff>critin .and. iter<maxiter)

        !$omp parallel do private(ix,ie,iy,id,jx,je,jy,jd,vcond,cf)
        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            ! compute the conditional expectaions
            vcond = 0.0d0
            do jx = 1,nx

                je = idmat(jx,1)
                jy = idmat(jx,2)
                jd = idmat(jx,3)
                vcond = vcond + Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*vmat0(:,jx)

            end do

            ! fit univeriate splines
            cf = spfit(invTk,vcond,rk,knotsk)

            !$omp parallel do private(ik,enow,ynow,beta,know,yterm,klow,khigh,kp,f0,df0,d2f0)
            do ik = 1,nk

                enow = Ge(ie)
                ynow = Gy(iy)
                beta = Gd(id)
                know = knotsk(ik)
                yterm = w0*ynow*enow + r0*know/THETA

                if (s1==0) then

                    klow = knotsk(1)
                    khigh = min(yterm,knotsk(nk))
                    if (nraflag) then
                        call nra(gmat0(ik,ix),klow,khigh,kp,cf,yterm,knotsk,beta)
                    else
                        call gss(klow,khigh,kp,cf,yterm,knotsk,beta)
                    end if
                    call vfuncsp(kp,cf,yterm,knotsk,beta,f0,df0,d2f0)
                    gmat1(ik,ix) = kp
                    vmat1(ik,ix) = -f0

                else

                    kp = gmat0(ik,ix)
                    call vfuncsp(kp,cf,yterm,knotsk,beta,f0,df0,d2f0)
                    gmat1(ik,ix) = kp
                    vmat1(ik,ix) = -f0

                end if

            end do ! loop for k
            !$omp end parallel do

        end do ! loop for x
        !$omp end parallel do

        diffg = maxval(maxval(abs(gmat1-gmat0),1),1)
        diffv = maxval(maxval(abs(vmat1-vmat0),1),1)
        diff  = diffv
        iter  = iter + 1

        ! diagnosis
        if (mod(iter,diagnum)==0) then

            write(*,"('  iteration ', I4, '  ||Tv-v|| = ', F10.5, '  ||Tg-g|| = ', F10.5)") iter, diffv, diffg

        end if

        if (diffg<1d-4 .and. s1==0) then
            s1 = 1
        elseif (s1>0 .and. s1<20) then
            s1 = s1+1
        elseif (s1>=20) then
            s1 = 0
        end if

        gmat0 = gmat1
        vmat0 = vmat1

    end do


    ! fit splines
    if (spliflag) then

        do ix = 1,nx

            cfmate(:,:,ix) = spfit(invTk,gmat0(:,ix),rk,knotsk)

        end do

    else

        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            ! compute the conditional expectaions
            vcond = 0.0d0
            do jx = 1,nx

                je = idmat(jx,1)
                jy = idmat(jx,2)
                jd = idmat(jx,3)
                vcond = vcond + Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*vmat0(:,jx)

            end do

            cfmate(:,:,ix) = spfit(invTk,vcond,rk,knotsk)

        end do

    end if


    call system_clock(ct1,cr)

    if (transmat==0) then
        ! iterarive method without transition matrix
        ! it is faster as we will not visit places with no population???
        diff = 1d+4
        iter = 0

        SOLVE = 0
        G = 0.0d0
        ! NOTE: need to initialize matrices as some points in the rectangle grid are never visited
        mpmat = 0.0d0
        kb = 0
        wb = 0.0d0

        do while(diff>critmu) ! .and. iter<maxiter)

            !$omp parallel do private(ix,ib,know,ie,iy,id,enow,ynow,beta,yterm,klow,khigh,kp,df0,d2f0)
            do ix = 1,nx

                do ib = 1,nb

                    ! if the distribution is greater than 0
                    if (mu0(ib,ix)>0.0d0) then

                        ! if the grid point is not visited before
                        if (SOLVE(ib,ix)==0) then

                            know = knotsb(ib)

                            if (spliflag) then

                                ! evaluate k' = g(k,e) by using splines
                                call speva(reshape(cfmate(:,:,ix),(/4,rk+1/)),know,rk,knotsk,kp,df0,d2f0)

                            else

                                ie = idmat(ix,1)
                                iy = idmat(ix,2)
                                id = idmat(ix,3)

                                enow = Ge(ie)
                                ynow = Gy(iy)
                                beta = Gd(id)
                                yterm = w0*ynow*enow + r0*know/THETA

                                klow = knotsk(1)
                                khigh = min(yterm,knotsk(nk))
                                if (nraflag) then
                                    call nra(know,klow,khigh,kp,cfmate(:,:,ix),yterm,knotsk,beta)
                                else
                                    call gss(klow,khigh,kp,cfmate(:,:,ix),yterm,knotsk,beta)
                                end if

                            end if

                            G(ib,ix) = kp
                            SOLVE(ib,ix) = 1

                        else

                            kp = G(ib,ix)

                        end if

                        mpmat(ib,ix) = kp*THETA
                        kb(ib,ix) = gridlookup2(kp,knotsb)
                        wb(ib,ix) = (knotsb(kb(ib,ix)+1)-kp)/(knotsb(kb(ib,ix)+1)-knotsb(kb(ib,ix)))

                    end if

                end do

            end do
            !$omp end parallel do

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

            diff = maxval(maxval(abs(mu1-mu0),1),1)
            iter = iter + 1
            ! diagnosis
            if (mod(iter,500)==0) then

                write(*,"('  iteration ', I6, '  ||mu1-mu0|| = ', F15.5, ' sum(mu1) = ', F15.5)") iter, diff, sum(sum(mu1,1),1)

            end if

            ! update distribution
            mu0 = mu1

        end do ! end mu loop

    else

        ! with transition matrix
        !$omp parallel do private(ix,ib,know,ie,iy,id,enow,ynow,beta,yterm,klow,khigh,kp,df0,d2f0)
        do ix = 1,nx

            do ib = 1,nb

                know = knotsb(ib)

                if (spliflag) then

                    ! evaluate k' = g(k,e) by using splines
                    call speva(reshape(cfmate(:,:,ix),(/4,rk+1/)),know,rk,knotsk,kp,df0,d2f0)

                else

                    ie = idmat(ix,1)
                    iy = idmat(ix,2)
                    id = idmat(ix,3)

                    enow = Ge(ie)
                    ynow = Gy(iy)
                    beta = Gd(id)
                    yterm = w0*ynow*enow + r0*know/THETA

                    klow = knotsk(1)
                    khigh = min(yterm,knotsk(nk))
                    if (nraflag) then
                        call nra(know,klow,khigh,kp,cfmate(:,:,ix),yterm,knotsk,beta)
                    else
                        call gss(klow,khigh,kp,cfmate(:,:,ix),yterm,knotsk,beta)
                    end if

                end if

                mpmat(ib,ix) = kp*THETA
                kb(ib,ix) = gridlookup2(kp,knotsb)
                wb(ib,ix) = (knotsb(kb(ib,ix)+1)-kp)/(knotsb(kb(ib,ix)+1)-knotsb(kb(ib,ix)))

            end do

        end do
        !$omp end parallel do

        ! allocate(aa(nb*nx,nb*nx))
        ! AA = 0.0d0
        ! print *, count(kb==1)
        allocate(csrColInd(nb*nx*nx*3-count(kb==1)),csrVal(nb*nx*nx*3-count(kb==1)))
        index = 0

        do ix = 1,nx

            ie = idmat(ix,1)
            iy = idmat(ix,2)
            id = idmat(ix,3)

            do ib = 1,nb

                csrRowPtr(nb*(ix-1)+ib) = index+1

                do jx = 1,nx

                    je = idmat(jx,1)
                    jy = idmat(jx,2)
                    jd = idmat(jx,3)

                    if (kb(ib,ix)==1) then

                        ! newborns & wb fraction
                        ! AA(nb*(ix-1)+ib,nb*(jx-1)+1) = wb(ib,ix)*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA + Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*(1.0d0-THETA)
                        index = index+1
                        csrVal(index) = wb(ib,ix)*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA + Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*(1.0d0-THETA)
                        csrColInd(index) = nb*(jx-1)+kb(ib,ix)

                    else
                        ! newborns
                        ! AA(nb*(ix-1)+ib,nb*(jx-1)+1) = Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*(1.0d0-THETA)
                        index = index+1
                        csrVal(index) = Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*(1.0d0-THETA)
                        csrColInd(index) = nb*(jx-1)+1

                        ! AA(nb*(ix-1)+ib,nb*(jx-1)+kb(ib,ix)) = wb(ib,ix)*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA
                        index = index+1
                        csrVal(index) = wb(ib,ix)*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA
                        csrColInd(index) = nb*(jx-1)+kb(ib,ix)

                    end if

                    ! AA(nb*(ix-1)+ib,nb*(jx-1)+kb(ib,ix)+1) = (1.0d0-wb(ib,ix))*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA
                    index = index+1
                    csrVal(index) = (1.0d0-wb(ib,ix))*Pe(ie,je)*Py(iy,jy)*Pd(id,jd)*THETA
                    csrColInd(index) = nb*(jx-1)+kb(ib,ix)+1

                end do

            end do

        end do

        csrRowPtr(nb*nx+1) = index+1
        !   Create CSR matrix
        ! i = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ZERO,nb*ne,nb*ne,csrRowPtr,csrRowPtr(2),csrColInd,csrVal)
        i = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,nb*nx,nb*nx,csrRowPtr,csrRowPtr(2),csrColInd,csrVal)
        !   Create matrix descriptor
        descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
        !   Analyze sparse matrix; chose proper kernels and workload balancing strategy
        info = MKL_SPARSE_OPTIMIZE(csrA)

        ! Option 1: iterative method
        if (transmat==1) then

            diff = 1d+4
            iter = 0

            do while(diff>critmu) ! .and. iter<maxiter)

                ! temp = matmul(transpose(AA),reshape(mu0,(/nb*nx/)))
                ! call dgemv('T', nb*nx, nb*nx, 1.0d0, AA, nb*nx, reshape(mu0,(/nb*nx/)), 1, 0.0d0, temp, 1)
                info = MKL_SPARSE_D_MV(SPARSE_OPERATION_TRANSPOSE,1.0d0,csrA,descrA,reshape(mu0,(/nb*nx/)),0.0d0,temp)

                mu1 = reshape(temp,(/nb,nx/))

                diff = maxval(maxval(abs(mu1-mu0),1),1)
                iter = iter + 1
                ! diagnosis
                if (mod(iter,500)==0) then

                    write(*,"('  iteration ', I6, '  ||mu1-mu0|| = ', F20.15, ' sum(mu1) = ', F20.15)") iter, diff, sum(sum(mu1,1),1)

                end if

                ! update distribution
                mu0 = mu1

            end do

        ! Option 2: lapack (dgeev)
        else if (transmat==2) then

            call eig(transpose(AA),mu0)

        ! Option 3: arpack (dnaupd and dneupd)
        else if (transmat==3) then

            call eigs(AA,csrA,descrA,1,'LM',critmu,mu0)

        end if

        !   Release internal representation of CSR matrix
        info = MKL_SPARSE_DESTROY(csrA)

    end if


    ! aggregate k
    m1 = 0.0d0
    do ix = 1,nx

        m1 = m1 + sum(knotsb*mu0(:,ix),1)

    end do

    print *, sum(mu0,1)
    call system_clock(ct2, cr)
    eptime = dble(ct2-ct1)/cr
    write(*,"('  Elasped time = ', F10.5)") eptime


end subroutine driverm


subroutine nra(xinit,xmin,xmax,x,cf,yterm,knotsk,beta)


    use mod_parameters, only: rk, critn
	real(8), intent (in) :: xinit, xmin, xmax, cf(:,:), yterm, knotsk(:), beta
	real(8), intent (out) :: x
    real(8) FL, DFL, FH, DFH, f0, df0, d2f0, NewtonStep, x0, xlow, xhigh, x1, diff
    integer iter


    call vfuncsp(xmin,cf,yterm,knotsk,beta,FL,DFL,D2F0)
    call vfuncsp(xmax,cf,yterm,knotsk,beta,FH,DFH,D2F0)

    if (DFL<0.0d0) then

        x0 = xmin

    elseif (DFH>0.0d0) then

        x0 = xmax

    else ! DFL>=0 amd DFH<=0

        x0    = xinit
        xlow  = xmin
        xhigh = xmax

        diff = 1d+4
        iter = 0

        do while (diff>critn)

            call vfuncsp(x0,cf,yterm,knotsk,beta,F0,DF0,D2F0)
            NewtonStep = df0/d2f0
            x1 = x0 - NewtonStep

            if (df0>0.0d0) then
                xlow = x0
            else
                xhigh = x0
            end if

            ! backtrack
            if (x1>xhigh) then
                x1 = (xlow+xhigh)/2.0d0
            elseif (x1<xlow) then
                x1 = (xlow+xhigh)/2.0d0
            end if

            diff = abs(x1-x0)
            iter = iter + 1
            ! diagnosis
            ! write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff
            x0 = x1

        end do

        ! check bounds
        if (x0<xmin) then
            x0 = xmin
        elseif (x0>xmax) then
            x0 = xmax
        end if

    end if

    x = x0


end subroutine nra


subroutine gss(xmin,xmax,x,cf,yterm,knotsk,beta)


    use mod_parameters, only: rk, critg
	real(8), intent (in) :: xmin, xmax, cf(:,:), yterm, knotsk(:), beta
	real(8), intent (out) :: x
	real(8) rg, a, b, c, d, z, fc, df0, d2f0, fd, diff
	integer iter


	rg = (3.0d0-sqrt(5.0d0))/2.0d0

	a = xmin
	b = xmax
	c = a + rg*(b-a)
    call vfuncsp(c,cf,yterm,knotsk,beta,fc,DF0,D2F0)
	d = a + (1.0d0-rg)*(b-a)
    call vfuncsp(d,cf,yterm,knotsk,beta,fd,DF0,D2F0)

!	crit = 1d-5
    diff = 1d+4
	iter = 0

	do while (diff>critg)

		if (fc>=fd) then

            z = c + (1.0d0-rg)*(b-c)
           	a = c
           	c = d
           	fc = fd
           	d = z
            call vfuncsp(d,cf,yterm,knotsk,beta,fd,DF0,D2F0)

        else

           	z = a + rg*(d-a)
           	b = d
           	d = c
           	fd = fc
           	c = z
            call vfuncsp(c,cf,yterm,knotsk,beta,fc,DF0,D2F0)

       	end if

        diff = d-c
        iter = iter+1

        ! diagnosis
        ! write(*,"('  iteration ', I4, '  New x = ', F10.5)") iter, diff

	end do

	x = c


end subroutine gss


subroutine vfuncsp(k0,cf,yterm,knotsk,beta,F0,DF0,D2F0)


    use mod_parameters, only: rk, SIGMA, THETA
	real(8), intent(in) :: k0, cf(:,:), yterm, knotsk(:), beta
    real(8), intent (out) :: F0, DF0, D2F0
	real(8) cnow, EV, EDV, ED2V


    cnow = yterm - k0
    call speva(cf,k0,rk,knotsk,EV,EDV,ED2V)

    if (SIGMA==1.0d0) then
        F0 = log(cnow) + BETA*THETA*ev
    else
        F0 = (cnow**(1.0d0-SIGMA)-1.0d0)/(1.0d0-SIGMA) + BETA*THETA*ev
    end if
    F0   = -F0 ! for gss
    DF0  = -cnow**(-SIGMA) + BETA*THETA*edv
    D2F0 = -SIGMA*cnow**(-SIGMA-1.0d0) + BETA*THETA*ed2v


end subroutine vfuncsp


subroutine calcdiststat(kgrid,mumat,thresholdwvec,sharewvec,gini)
    ! wealth distribution
	! calculating threshold
    real(8), intent(in) :: kgrid(:), mumat(:,:)
    integer, intent(out) :: thresholdwvec(:)
    real(8), intent(out) :: sharewvec(:), gini
    integer ne, nk, ie, ik
    real(8), allocatable :: muk(:)
	real(8) cummmk, wgtik(7), negW, ShareW5(8)


	nk = size(kgrid,1)
    ne = size(mumat,2)

    allocate(muk(nk))

    muk = sum(mumat,2)

    call ginicoef(muk,kgrid,nk,gini) !,ShareW5,negW)

	thresholdwvec = 0         ! index
	cummmk = 0.0d0
	wgtik = 0.0d0           ! calculating division rule for threshold

	do ik = 1,nk
		cummmk = cummmk + muk(ik)
		if (cummmk >= 0.2d0 .and. thresholdwvec(1) == 0 ) then
			thresholdwvec(1) = ik
			wgtik(1) = 1.0d0 - (cummmk - 0.2d0)/(muk(ik))
		elseif (cummmk >= 0.4d0 .and. thresholdwvec(2) == 0 ) then
			thresholdwvec(2) = ik
			wgtik(2) = 1.0d0 - (cummmk - 0.4d0)/(muk(ik))
		elseif (cummmk >= 0.6d0 .and. thresholdwvec(3) == 0 ) then
			thresholdwvec(3) = ik
			wgtik(3) = 1.0d0 - (cummmk - 0.6d0)/(muk(ik))
		elseif (cummmk >= 0.8d0 .and. thresholdwvec(4) == 0 ) then
			thresholdwvec(4) = ik
			wgtik(4) = 1.0d0 - (cummmk - 0.8d0)/(muk(ik))
        elseif (cummmk >= 0.9d0 .and. thresholdwvec(5) == 0 ) then
			thresholdwvec(5) = ik
			wgtik(5) = 1.0d0 - (cummmk - 0.9d0)/(muk(ik))
        elseif (cummmk >= 0.95d0 .and. thresholdwvec(6) == 0 ) then
			thresholdwvec(6) = ik
			wgtik(6) = 1.0d0 - (cummmk - 0.95d0)/(muk(ik))
        elseif (cummmk >= 0.99d0 .and. thresholdwvec(7) == 0 ) then
			thresholdwvec(7) = ik
			wgtik(7) = 1.0d0 - (cummmk - 0.99d0)/(muk(ik))
		end if
	end do

	ShareW5 = 0.0d0

	do ik = 1,nk
		do ie = 1,ne
			if (ik < thresholdwvec(1)) then
				ShareW5(1) = ShareW5(1) + kgrid(ik)*mumat(ik,ie)
			elseif (ik == thresholdwvec(1)) then
				ShareW5(1) = ShareW5(1) + kgrid(ik)*mumat(ik,ie)*wgtik(1)
				ShareW5(2) = ShareW5(2) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(1))
			elseif (ik > thresholdwvec(1) .and. ik < thresholdwvec(2)) then
				ShareW5(2) = ShareW5(2) + kgrid(ik)*mumat(ik,ie)

			elseif (ik == thresholdwvec(2)) then
				ShareW5(2) = ShareW5(2) + kgrid(ik)*mumat(ik,ie)*wgtik(2)
				ShareW5(3) = ShareW5(3) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(2))
			elseif (ik > thresholdwvec(2) .and. ik < thresholdwvec(3)) then
				ShareW5(3) = ShareW5(3) + kgrid(ik)*mumat(ik,ie)

			elseif (ik == thresholdwvec(3)) then
				ShareW5(3) = ShareW5(3) + kgrid(ik)*mumat(ik,ie)*wgtik(3)
				ShareW5(4) = ShareW5(4) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(3))
			elseif (ik > thresholdwvec(3) .and. ik < thresholdwvec(4)) then
				ShareW5(4) = ShareW5(4) + kgrid(ik)*mumat(ik,ie)

			elseif (ik == thresholdwvec(4)) then
				ShareW5(4) = ShareW5(4) + kgrid(ik)*mumat(ik,ie)*wgtik(4)
				ShareW5(5) = ShareW5(5) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(4))
            elseif (ik > thresholdwvec(4) .and. ik < thresholdwvec(5)) then
				ShareW5(5) = ShareW5(5) + kgrid(ik)*mumat(ik,ie)

			elseif (ik == thresholdwvec(5)) then
				ShareW5(5) = ShareW5(5) + kgrid(ik)*mumat(ik,ie)*wgtik(5)
				ShareW5(6) = ShareW5(6) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(5))
            elseif (ik > thresholdwvec(5) .and. ik < thresholdwvec(6)) then
				ShareW5(6) = ShareW5(6) + kgrid(ik)*mumat(ik,ie)

			elseif (ik == thresholdwvec(6)) then
				ShareW5(6) = ShareW5(6) + kgrid(ik)*mumat(ik,ie)*wgtik(6)
				ShareW5(7) = ShareW5(7) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(6))
            elseif (ik > thresholdwvec(6) .and. ik < thresholdwvec(7)) then
				ShareW5(7) = ShareW5(7) + kgrid(ik)*mumat(ik,ie)

			elseif (ik == thresholdwvec(7)) then
				ShareW5(7) = ShareW5(7) + kgrid(ik)*mumat(ik,ie)*wgtik(7)
				ShareW5(8) = ShareW5(8) + kgrid(ik)*mumat(ik,ie)*(1.0d0-wgtik(7))
			else
				ShareW5(8) = ShareW5(8) + kgrid(ik)*mumat(ik,ie)

			end if
		end do
	end do

	sharewvec(1) = 100.0d0*ShareW5(1)/sum(ShareW5)
	sharewvec(2) = 100.0d0*ShareW5(2)/sum(ShareW5)
	sharewvec(3) = 100.0d0*ShareW5(3)/sum(ShareW5)
	sharewvec(4) = 100.0d0*ShareW5(4)/sum(ShareW5)
	sharewvec(5) = 100.0d0*sum(ShareW5(5:8))/sum(ShareW5)
    sharewvec(6) = 100.0d0*ShareW5(6)/sum(ShareW5)
    sharewvec(7) = 100.0d0*ShareW5(7)/sum(ShareW5)
    sharewvec(8) = 100.0d0*ShareW5(8)/sum(ShareW5)


end subroutine calcdiststat


subroutine eig(AA,mu0)


    real(8), intent(in) :: AA(:,:)
    real(8), intent(out) :: mu0(:,:)
    integer info, lda, ldvr, lwork, n, nb, nx
    integer, parameter :: nblock = 2048
    ! .. Local Arrays ..
    real(8), allocatable :: vr(:,:), wi(:), work(:), wr(:)
    real(8) dummy(1,1)


    nb = size(mu0,1)
    nx = size(mu0,2)
    n = nb*nx
    lda = n
    ldvr = n
    allocate(vr(ldvr,n),wi(n),wr(n))

    lwork = -1
    call dgeev('N', 'V', n, AA, lda, wr, wi, dummy, 1, vr, ldvr, dummy, lwork, info)

    ! Make sure that there is enough workspace for block size nb.
    lwork = max((nblock+2)*n,nint(dummy(1,1)))
    print *, lwork, nint(dummy(1,1))
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

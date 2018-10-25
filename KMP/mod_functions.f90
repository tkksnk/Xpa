module mod_functions


implicit none


contains


! subroutine ginicoef(dist,grid,ngrid,ns,gini) !,share,negW)
subroutine ginicoef(dist,grid,ngrid,gini) !,share,negW)

    ! Calculate Gini and Lorenz curve

    ! Inputs
    ! dist: measure
    ! grid: grid
    ! ngrid: number of grid points
    ! ns: sharenum

    ! Output
    ! gini: Gini Index
    ! share: Share in bin
    ! negW: fraction with non-positive
    implicit none

    integer, intent(in) :: ngrid !, ns
    real(8), intent(in) :: dist(:), grid(:)
    real(8), intent(out) :: gini !, negW
    ! real(8), intent(out), dimension(ns) :: share

    integer :: i1, j1
    real(8) :: aaa, bbb !, ccc, ddd
    real(8), dimension(ngrid) :: xaxis, a1, yaxis
    real(8), dimension(ngrid-1) :: dc, db
    ! real(8), dimension(ns) :: sharecum

    ! ddd = ns

    ! Cacluate the cdf
    xaxis=0.0d0
    do i1 = 1,ngrid
        xaxis(i1)=sum(dist(1:i1))
    end do

    ! Generate the Lorenz curve (yaxis)
    a1=0.0d0
    yaxis=0.0d0
    a1=grid*dist
    aaa=sum(a1)
    do i1 = 1,ngrid
        bbb=sum(a1(1:i1))
        yaxis(i1)=bbb/aaa
    end do
    ! Calcualte the gini coefficient of wealth
    do i1 = 1,ngrid-1
        dc(i1)=xaxis(i1+1)-xaxis(i1)
        db(i1)=(yaxis(i1+1)+yaxis(i1))*0.5d0
    end do

    bbb=sum(dc*db)
    aaa=0.5d0-bbb

    gini=2.0d0*aaa

    ! sharecum=0.0d0
    ! do j1 = 1,ns
    !     ccc = j1/ddd
    !     sharecum(j1) = lini(xaxis,yaxis,ccc)
    ! end do
    !
    ! share = 0.0d0
    ! share(1) = sharecum(1)
    ! do j1 = 1,ns-2
    !     share(j1+1) = sharecum(j1+1) - sharecum(j1)
    ! end do
    ! share(ns) = 1 - sharecum(ns-1)
    !
    ! ! share of nonpositive holders
    ! negW = 0.0d0
    ! do i1 = 1,ngrid
    !     if (grid(i1) <= 0) then
    !         negW = negW + dist(i1)
    !     end if
    ! end do

end subroutine ginicoef


subroutine Beta2mat(BetaK,knotsm,mpmat)


    real(8), intent(in) :: BetaK(:,:), knotsm(:)
    real(8), intent(out) :: mpmat(:,:)
    real(8) mnow
    integer nz, nm, iz, im


    nm = size(mpmat,1)
    nz = size(mpmat,2)

    do iz = 1,nz

        do im = 1,nm

            mnow = knotsm(im)
            mpmat(im,iz) = exp( BetaK(iz,1) + BetaK(iz,2)*log(mnow) )

        end do

    end do


end subroutine Beta2mat


subroutine calcPeKMP(pgg,pbb,pgg00,pgg10,pgb00,pgb10,pbg00,pbg10,pbb00,pbb10,rhoy,sigy,dbar,deps,Pez,Gy,Py,Gd,Pd)


    real(8), intent(in) :: pgg, pbb, pgg00, pgg10, pgb00, pgb10, pbg00, pbg10, pbb00, pbb10, rhoy, sigy, dbar, deps
    real(8), intent(out) :: Pez(:,:,:), Gy(:), Py(:,:), Gd(:), Pd(:,:)
    real(8) Pegg(2,2), Pegb(2,2), Pebg(2,2), Pebb(2,2)
    integer ny, nd, id

    ny = size(Gy,1)
    nd = size(Gd,1)

    Pegg(1,:) = (/pgg00, 1.0d0-pgg00/)
    Pegg(2,:) = (/pgg10, 1.0d0-pgg10/)
    Pegb(1,:) = (/pgb00, 1.0d0-pgb00/)
    Pegb(2,:) = (/pgb10, 1.0d0-pgb10/)
    Pebg(1,:) = (/pbg00, 1.0d0-pbg00/)
    Pebg(2,:) = (/pbg10, 1.0d0-pbg10/)
    Pebb(1,:) = (/pbb00, 1.0d0-pbb00/)
    Pebb(2,:) = (/pbb10, 1.0d0-pbb10/)

    Pegg = Pegg*pgg
    Pegb = Pegb*(1.0d0-pgg)
    Pebg = Pebg*(1.0d0-pbb)
    Pebb = Pebb*pbb

    Pez(:,:,1) = Pegg
    Pez(:,:,2) = Pegb
    Pez(:,:,3) = Pebg
    Pez(:,:,4) = Pebb

    if (ny==1) then

        Gy = 1.0d0
        Py = 1.0d0

    else

        call rouwenhorst(rhoy,sigy,ny,Gy,Py)
        Gy = exp(Gy)

    end if

    if (nd==1) then

        Gd = 0.99d0 !89975d0
        Pd = 1.0d0

    else

        Gd = linspace(dbar-deps,dbar+deps,nd)
        Pd = 0.0d0
        do id = 1,nd

            Pd(id,id) = 1.0d0

        end do

    end if


end subroutine calcPeKMP


subroutine calcPeKS(pgg,pbb,ug,ub,DurationUg,DurationUb,corr,Pez)


    real(8), intent(in) :: pgg, pbb, ug, ub, DurationUg, DurationUb, corr
    real(8), intent(out) :: Pez(:,:,:)
    real(8) pgg00, pgg10, pbb00, pbb10, pgb00, pgb10, pbg00, pbg10, Pegg(2,2), Pegb(2,2), Pebg(2,2), Pebb(2,2)

    pgg00 = 1.0d0-1.0d0/DurationUg
    pgg10 = 1.0d0/(1.0d0-ug)*(ug-ug*pgg00)

    pbb00 = 1.0d0-1.0d0/DurationUb
    pbb10 = 1.0d0/(1.0d0-ub)*(ub-ub*pbb00)

    pgb00 = (1.0d0+corr)*pbb00
    pgb10 = 1.0d0/(1.0d0-ug)*(ub-ug*pgb00)

    pbg00 = (1.0d0-corr)*pgg00
    pbg10 = 1.0d0/(1.0d0-ub)*(ug-ub*pbg00)

    Pegg(1,:) = (/pgg00, 1.0d0-pgg00/)
    Pegg(2,:) = (/pgg10, 1.0d0-pgg10/)
    Pegb(1,:) = (/pgb00, 1.0d0-pgb00/)
    Pegb(2,:) = (/pgb10, 1.0d0-pgb10/)
    Pebg(1,:) = (/pbg00, 1.0d0-pbg00/)
    Pebg(2,:) = (/pbg10, 1.0d0-pbg10/)
    Pebb(1,:) = (/pbb00, 1.0d0-pbb00/)
    Pebb(2,:) = (/pbb10, 1.0d0-pbb10/)

    Pegg = Pegg*pgg
    Pegb = Pegb*(1.0d0-pgg)
    Pebg = Pebg*(1.0d0-pbb)
    Pebb = Pebb*pbb

    Pez(:,:,1) = Pegg
    Pez(:,:,2) = Pegb
    Pez(:,:,3) = Pebg
    Pez(:,:,4) = Pebb


end subroutine calcPeKS


function kron(amat,bmat) result(cmat)

    real(8), intent(in) :: amat(:,:), bmat(:,:)
    integer na, ma, nb, mb, nc, mc, ia, ja, ib, jb
    real(8), allocatable :: cmat(:,:)
    real(8) a0, b0

    na = size(amat,1)
    ma = size(amat,2)
    nb = size(bmat,1)
    mb = size(bmat,2)
    nc = na*nb
    mc = ma*mb
    allocate(cmat(nc,mc))

    do ia = 1,na

        do ja = 1,ma

            a0 = amat(ia,ja)

            do ib = 1,nb

                do jb = 1,mb

                    b0 = bmat(ib,jb)
                    cmat(nb*(ia-1)+ib,mb*(ja-1)+jb) = a0*b0

                end do

            end do

        end do

    end do

end function kron


function calcmu(PIe) result(mue)


    ! use mod_parameters
    ! use mod_calcss, only: eig
    real(8), intent(in) :: PIe(:,:)
    real(8), allocatable :: Pemat(:,:), mue(:)
    integer ie, ne, i
    integer info, lda, ldvr, lwork
    integer, parameter :: nblock = 2048
    ! .. Local Arrays ..
    real(8), allocatable :: vr(:,:), wi(:), work(:), wr(:)
    real(8) dummy(1,1)

    ne = size(PIe,1)
    allocate(Pemat(ne,ne),mue(ne))

    lda = ne
    ldvr = ne
    allocate(vr(ldvr,ne),wi(ne),wr(ne))
    lwork = (nblock+2)*ne !max((nblock+2)*n,nint(dummy(1,1)))
    ! print *, lwork, nint(dummy(1,1))
    allocate(work(lwork))

    ! stationary distribution of e
    ! Compute the eigenvalues and right eigenvectors of A
    call dgeev('N', 'V', ne, transpose(PIe), lda, wr, wi, dummy, 1, vr, ldvr, work, lwork, info)
    mue = reshape(vr(:,maxloc(wr))/sum(vr(:,maxloc(wr))),(/ne/))

    ! call eig(transpose(PIe),mue)
    ! do ie = 1,ne
    !     Pemat(ie,ie) = 1.0d0
    ! end do
    !
    ! do i = 1,10000
    !     Pemat = matmul(Pemat,PIe)
    ! end do
    !
    ! mue = reshape(Pemat(1,:), (/ne/))
    ! mue = mue/sum(mue)


end function calcmu


subroutine qnwnorm(n,mu,sigma,x,w)


    integer, intent(in) :: n
    real(8), intent(in) :: mu, sigma
    real(8), intent(out) :: x(n), w(n)
    real(8) pi, pim4, z, p1, p2, p3, pp, z1
    integer maxit, m, its, i, j

    pi = 3.141592653589793d0
    maxit = 100
    pim4 = 1.0d0/pi**0.25d0
    m = max(floor(dble((n+1)/2.0d0)),0)

    do i = 1,m

        ! Reasonable starting values
        if (i==1) then
            z = sqrt(dble(2*n+1))-1.85575d0*(dble(2*n+1)**(-1.0d0/6.0d0))
        elseif (i==2) then
            z = z-1.14d0*(dble(n)**0.426d0)/z
        elseif (i==3) then
            z = 1.86d0*z+0.86d0*x(1)
        elseif (i==4) then
            z = 1.91d0*z+0.91d0*x(2)
        else
            z = 2.0d0*z+x(i-2)
        end if

        ! root finding iterations
        its = 0

        do while (its<maxit)

            its = its+1
            p1 = pim4
            p2 = 0.0d0

            do j = 1,n

                p3 = p2
                p2 = p1
                p1 = z*sqrt(2.0d0/dble(j))*p2-sqrt(dble(j-1)/dble(j))*p3

            end do

            pp = sqrt(dble(2*n))*p2
            z1 = z
            z  = z1-p1/pp

            if (abs(z-z1)<1d-14) then
                exit
            end if

        end do

        !if its>=maxit
        !   error('failure to converge in qnwnorm1')
        !end

        x(n+1-i) = z
        x(i) = -z
        w(i) = 2.0d0/(pp*pp)
        w(n+1-i) = w(i)

    end do

    w = w/sqrt(pi)
    x = x*sqrt(2.0d0)

    x = mu + sigma*x


end subroutine qnwnorm


subroutine rouwenhorst(rho,sigmas,znum,Z,PI)
! Creates a discrete approximation to a first order autoregressive process
! with serial correlation p + q - 1.  See Rouwenhorst (1995) (pps. 325-329)
! in Cooley <i>Frontiers of Business Cycle Research</i> Princeton.


    ! declare variables
    integer, intent(in) :: znum
    real(8), intent(in) :: rho, sigmas
    real(8), intent(out) :: Z(znum)
    real(8), intent(out) :: PI(znum,znum)

    real(8) p, q, zvar, epsilon
    real(8), allocatable :: h1(:,:), h2(:,:), h3(:,:), h4(:,:), h(:,:), hlag(:,:)
    integer i


    p = (rho + 1.0d0)/2.0d0
    q = p

    allocate(hlag(1,1))
    hlag(1,1) = 1.0d0

    do i = 2,znum

	    !    vec0 = zeros(i-1,1);
        allocate(h1(i,i), h2(i,i), h3(i,i), h4(i,i), h(i,i))
        h1 = 0.0d0
        h2 = 0.0d0
        h3 = 0.0d0
        h4 = 0.0d0
        h1(1:i-1,1:i-1) = hlag
        h2(1:i-1,2:i)   = hlag
        h3(2:i,1:i-1)   = hlag
        h4(2:i,2:i)     = hlag
        h = p*h1 + (1-p)*h2 + (1-q)*h3 + q*h4
        h(2:i-1,:) = h(2:i-1,:)/2

        deallocate(hlag)
        allocate(hlag(i,i))
        hlag = h

        deallocate(h1, h2, h3, h4, h)
	    !    h = p*[hlag vec0; vec0' 0] + (1-p)*[vec0 hlag; 0 vec0'] + ...
	    ! 	    (1- q)*[vec0' 0; hlag vec0] + q*[0 vec0'; vec0 hlag]
	    !    h(2:i-1,:) = h(2:i-1,:)./2
        !    hlag = h

    end do

    PI = hlag
!    PI = 1.0d0

    ! symmetrically and evenly spaced between [-epsilon, epsilon] with h elements.
    ! When p = q, then then variance of shock is epsilon^2/(h-1).
    zvar = (sigmas**2)/(1.0d0 - rho**2)
    epsilon = sqrt((znum - 1.0d0)*zvar)

    Z = linspace(-epsilon, epsilon, znum)


end subroutine rouwenhorst


function linspace(a,b,m) result(x)


    integer, intent(in) :: m
    real(8), intent(in) :: a, b
    real(8), dimension(m) :: x
    integer i

    x(1) = a

    do i = 2,m

        x(i) = x(i-1) + (b-a)/(m-1)

    end do


end function linspace


function logspace(a,b,m) result(x)


    integer, intent(in) :: m
    real(8), intent(in) :: a, b
    real(8), dimension(m) :: x
    integer i

    x = linspace(log(10.0d0**a),log(10.0d0**b),m)
    x = exp(x)


end function logspace


! NOTE: 070718 from mod_utils.f90 in the JSY code
function gridlookup2(x0,xgrid) result(ix)


    real(8), intent(in) :: x0, xgrid(:)
    integer ix, jx, nx


	nx = size(xgrid,1)

    ! NOTE: modified on 030918
    ! ix = 0
    ! do jx = 1,nx
    !     if (xgrid(jx)>=x0) exit
    !     ix = ix+1
    ! end do
    ix = 0
    do jx = 1,nx
        ! if (xgrid(jx)>=x0) exit
        if (x0<=xgrid(jx)) exit
        ix = ix+1
    end do

    ix = min(max(1,ix),nx-1)


end function gridlookup2


function gridlookup(nb,mbgrid,b0) result(k)


    integer, intent(in) :: nb
    real(8), intent(in) :: mbgrid(nb), b0
    integer i, k

    k = 0
    do i = 1,nb
        if (mbgrid(i)>b0) exit
        k = k+1
    end do
!    k = sum(mbgrid<=b0)
    k = max(k,1)
    k = min(k,nb-1)


end function gridlookup


function cumsum(A) result(cumA)


    real(8), intent(in) :: A(:,:)
    real(8), allocatable :: cumA(:,:)
    integer i,j,n

    n = size(A,1)

    allocate(cumA(n,n))

    cumA = 0.0d0

    do i = 1,n

        do j = 1,n

            cumA(i,j) = sum(A(i,1:j),1)

        end do

    end do


end function cumsum


end module mod_functions

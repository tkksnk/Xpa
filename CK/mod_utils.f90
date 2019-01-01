module mod_utils
! NOTE: 030308 soubroutines readwritekappa and readwriteVmat are added

implicit none


contains


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


! function genzvec(simTT,RHOZ,SDINOVZ) result(zvec)
!
!
!     integer, intent(in) :: simTT
!     real(8), intent(in) :: RHOZ, SDINOVZ !pzij(:,:)
!     ! real(8), allocatable :: cumpz(:,:)
!     real(8) randn
!     integer seedsize, time
!     integer, allocatable :: seed(:), zvec(:)
!
!
!     ! nz = size(pzij,1)
!     ! allocate(cumpz(nz,nz))
!
!     call random_seed(size=seedsize)
!     allocate(seed(seedsize),zvec(simTT))
!     call random_seed(get=seed)
!
!     print *, "Size of seed array is", seedsize
!     call random_seed(put=seed)
!
!     ! cumpz = cumsum(pzij)
!     ! izvec(1) = ceiling(dble(nz)/2.0d0)
!     zvec(1) = 1.0d0
!     do time = 1,simTT-1
!
!         znow = zvec(time)
!         call random_number(randn) ! normal?
!         zvec(time+1) = exp(RHOZ*log(znow) + SDINOVZ*randn)
!
!     end do
!
!
! end function genzvec


function genizvec(simTT,pzij) result(izvec)


    integer, intent(in) :: simTT
    real(8), intent(in) :: pzij(:,:)
    real(8), allocatable :: cumpz(:,:)
    real(8) rand
    integer seedsize, time, nz
    integer, allocatable :: seed(:), izvec(:)


    nz = size(pzij,1)
    allocate(cumpz(nz,nz))

    call random_seed(size=seedsize)
    allocate(seed(seedsize),izvec(simTT))
    call random_seed(get=seed)

    print *, "Size of seed array is", seedsize
    call random_seed(put=seed)

    cumpz = cumsum(pzij)
    izvec(1) = ceiling(dble(nz)/2.0d0)
    do time = 1,simTT-1

        call random_number(rand)
        izvec(time+1) = count(rand-cumpz(izvec(time),:)>=0)
        izvec(time+1) = min(izvec(time+1)+1,nz)

    end do


end function genizvec


subroutine readwritekappa(kappakp,kappaw,kappar,kappal,ckappakp,ckappaw,ckappar,ckappal,conzflag,readflag)


    real(8), intent(out) :: kappakp(:,:), kappaw(:,:), kappar(:,:), kappal(:,:)
    integer, intent(in) :: conzflag, readflag
    character(40), intent(in) :: ckappakp, ckappaw, ckappar, ckappal
    integer iz, nz


    if (conzflag==1) then

        ! 23 May 2018: kappakp shape has been altered.
        !nz = size(kappakp,1)
		nz = size(kappakp,2)

        do iz = 1,nz
            if (readflag==1) then
                open(410, file=ckappakp)
                open(420, file=ckappaw)
                open(430, file=ckappar)
                open(440, file=ckappal)
                ! 23 May 2018: kappakp shape has been altered.
                !read(410,*) kappakp(iz,1)
                !read(410,*) kappakp(iz,2)
                !read(420,*) kappaw(iz,1)
                !read(420,*) kappaw(iz,2)
                !read(430,*) kappar(iz,1)
                !read(430,*) kappar(iz,2)
                !read(440,*) kappal(iz,1)
                !read(440,*) kappal(iz,2)
                read(410,*) kappakp(1,iz)
                read(410,*) kappakp(2,iz)
                read(420,*) kappaw(1,iz)
                read(420,*) kappaw(2,iz)
                read(430,*) kappar(1,iz)
                read(430,*) kappar(2,iz)
                read(440,*) kappal(1,iz)
                read(440,*) kappal(2,iz)
            else
                open(410, file="./kappakp.txt")
                open(420, file="./kappaw.txt")
                open(430, file="./kappar.txt")
                open(440, file="./kappal.txt")
                ! 23 May 2018: kappakp shape has been altered.
                !write(410,*) kappakp(iz,1)
                !write(410,*) kappakp(iz,2)
                !write(420,*) kappaw(iz,1)
                !write(420,*) kappaw(iz,2)
                !write(430,*) kappar(iz,1)
                !write(430,*) kappar(iz,2)
                !write(440,*) kappal(iz,1)
                !write(440,*) kappal(iz,2)
                write(410,*) kappakp(1,iz)
                write(410,*) kappakp(2,iz)
                write(420,*) kappaw(1,iz)
                write(420,*) kappaw(2,iz)
                write(430,*) kappar(1,iz)
                write(430,*) kappar(2,iz)
                write(440,*) kappal(1,iz)
                write(440,*) kappal(2,iz)
            end if
        end do

    elseif (conzflag==3) then

        if (readflag==1) then
            open(410, file=ckappakp)
            open(420, file=ckappaw)
            open(430, file=ckappar)
            open(440, file=ckappal)
            read(410,*) kappakp(1,1)
            read(410,*) kappakp(2,1)
            read(410,*) kappakp(3,1)
            read(410,*) kappakp(4,1)
            read(420,*) kappaw(1,1)
            read(420,*) kappaw(2,1)
            read(420,*) kappaw(3,1)
            read(420,*) kappaw(4,1)
            read(430,*) kappar(1,1)
            read(430,*) kappar(2,1)
            read(430,*) kappar(3,1)
            read(430,*) kappar(4,1)
            read(440,*) kappal(1,1)
            read(440,*) kappal(2,1)
            read(440,*) kappal(3,1)
            read(440,*) kappal(4,1)
        else
            open(410, file="./kappakp.txt")
            open(420, file="./kappaw.txt")
            open(430, file="./kappar.txt")
            open(440, file="./kappal.txt")
            write(410,*) kappakp(1,1)
            write(410,*) kappakp(2,1)
            write(410,*) kappakp(3,1)
            write(410,*) kappakp(4,1)
            write(420,*) kappaw(1,1)
            write(420,*) kappaw(2,1)
            write(420,*) kappaw(3,1)
            write(420,*) kappaw(4,1)
            write(430,*) kappar(1,1)
            write(430,*) kappar(2,1)
            write(430,*) kappar(3,1)
            write(430,*) kappar(4,1)
            write(440,*) kappal(1,1)
            write(440,*) kappal(2,1)
            write(440,*) kappal(3,1)
            write(440,*) kappal(4,1)
        end if

    else
        if (readflag==1) then
            open(410, file=ckappakp)
            open(420, file=ckappaw)
            open(430, file=ckappar)
            open(440, file=ckappal)
            read(410,*) kappakp(1,1)
            read(410,*) kappakp(2,1)
            read(410,*) kappakp(3,1)
            read(420,*) kappaw(1,1)
            read(420,*) kappaw(2,1)
            read(420,*) kappaw(3,1)
            read(430,*) kappar(1,1)
            read(430,*) kappar(2,1)
            read(430,*) kappar(3,1)
            read(440,*) kappal(1,1)
            read(440,*) kappal(2,1)
            read(440,*) kappal(3,1)
        else
            open(410, file="./kappakp.txt")
            open(420, file="./kappaw.txt")
            open(430, file="./kappar.txt")
            open(440, file="./kappal.txt")
            write(410,*) kappakp(1,1)
            write(410,*) kappakp(2,1)
            write(410,*) kappakp(3,1)
            write(420,*) kappaw(1,1)
            write(420,*) kappaw(2,1)
            write(420,*) kappaw(3,1)
            write(430,*) kappar(1,1)
            write(430,*) kappar(2,1)
            write(430,*) kappar(3,1)
            write(440,*) kappal(1,1)
            write(440,*) kappal(2,1)
            write(440,*) kappal(3,1)
        end if

    end if

    close(410)
    close(420)
    close(430)
    close(440)


end subroutine readwritekappa


! subroutine readwritekappa0(kappakp0,kappaw0,kappar0,kappal0,ckappakp0,ckappaw0,ckappar0,ckappal0,readflag)
!
!
!     real(8), intent(out) :: kappakp0(:), kappaw0(:), kappar0(:), kappal0(:)
!     integer, intent(in) :: readflag
!     character(40), intent(in) :: ckappakp0, ckappaw0, ckappar0, ckappal0
!     integer iz, nz
!
!
!     if (readflag==1) then
!         open(410, file=ckappakp0)
!         open(420, file=ckappaw0)
!         open(430, file=ckappar0)
!         open(440, file=ckappal0)
!         read(410,*) kappakp0(1)
!         read(410,*) kappakp0(2)
!         read(410,*) kappakp0(3)
!         read(420,*) kappaw0(1)
!         read(420,*) kappaw0(2)
!         read(420,*) kappaw0(3)
!         read(430,*) kappar0(1)
!         read(430,*) kappar0(2)
!         read(430,*) kappar0(3)
!         read(440,*) kappal0(1)
!         read(440,*) kappal0(2)
!         read(440,*) kappal0(3)
!     else
!         open(410, file="./kappakp0.txt")
!         open(420, file="./kappaw0.txt")
!         open(430, file="./kappar0.txt")
!         open(440, file="./kappal0.txt")
!         write(410,*) kappakp0(1)
!         write(410,*) kappakp0(2)
!         write(410,*) kappakp0(3)
!         write(420,*) kappaw0(1)
!         write(420,*) kappaw0(2)
!         write(420,*) kappaw0(3)
!         write(430,*) kappar0(1)
!         write(430,*) kappar0(2)
!         write(430,*) kappar0(3)
!         write(440,*) kappal0(1)
!         write(440,*) kappal0(2)
!         write(440,*) kappal0(3)
!     end if
!
!     close(410)
!     close(420)
!     close(430)
!     close(440)
!
!
! end subroutine readwritekappa0


subroutine readwriteVmat(Vmat,readflag)


    real(8), intent(out) :: Vmat(:,:,:,:)
    integer, intent(in) :: readflag
    integer iz, ie, im, ia, nz, ne, nm, na


    na = size(Vmat,1)
    nm = size(Vmat,2)
    ne = size(Vmat,3)
    nz = size(Vmat,4)

    open(500, file="Vmat.txt")
    do iz = 1,nz
        do ie = 1,ne
            do im = 1,nm
                do ia = 1,na

                    if (readflag==1) then
                        read(500,*) Vmat(ia,im,ie,iz)
                    else
                        write(500,*) Vmat(ia,im,ie,iz)
                    end if

                end do
            end do
        end do
    end do
    close(500)


end subroutine readwriteVmat


subroutine findiw(x,xx,ix,wx)
! NOTE: For linear interpolation. How is it different from findjw?
    real(8), intent(in) :: x, xx(:)
    integer, intent(out) :: ix
    real(8), intent(out) :: wx
    integer nx


    nx = size(xx,1)

    if (x<=xx(1)) then
        ix = 1
    elseif (x>=xx(nx)) then
        ix = nx-1
    else
        do ix=1,nx-1

            if (x<xx(ix+1)) exit

        end do
        ! jx = 2
        ! do while (jx<=nx)
        !
        !     if (x<xx(jx)) then
        !         ix=jx-1
        !         jx=nx
        !     end if
        !
        ! jx=jx+1
        ! end do
    end if

    wx = (xx(ix+1)-x)/(xx(ix+1)-xx(ix))


end subroutine findiw


subroutine findjw(ap,kgrid,jk,wk)
! NOTE: For linear interpolation in distribution. How is it different from findiw?
    real(8), intent(in) :: ap, kgrid(:)
    integer, intent(out) :: jk
    real(8), intent(out) :: wk
    integer kk, nk


    nk = size(kgrid,1)

    if (ap<=kgrid(1)) then
        jk = 1
        wk = 1.0d0
    elseif (ap>=kgrid(nk)) then
        jk = nk-1
        wk = 0.0d0
    else
        jk = 0
        do kk = 1,nk
            ! if (kgrid(kk)>=ap) exit
            if (ap<=kgrid(kk)) exit
            jk = jk+1
        end do
        ! jk = min(max(1,jk),nk-1)
        wk = (kgrid(jk+1)-ap)/(kgrid(jk+1)-kgrid(jk))
    end if


end subroutine findjw


subroutine mpicounter(nproc,rank,ngrid,ncell,mystart,myend,mycount,acount,displacement)

    ! NOTE: ncell is the number of cells calculated by each iteration iv
	integer, intent(in) :: nproc, rank, ngrid, ncell
	integer, intent(out) :: mystart, myend, mycount, acount(:), displacement(:)
	! local variables for parallel code
    integer remainder, nlocalstates, j
    integer countvec(nproc), startvec(nproc), endvec(nproc)


	remainder = modulo(ngrid,nproc)
	nlocalstates = ngrid/nproc
	countvec = nlocalstates ! for 1,...,nproc

	do j = 1,remainder ! less than nproc
		countvec(j) = countvec(j) + 1
	end do
	startvec(1) = 1
	do j = 2,nproc
		startvec(j) = startvec(j-1) + countvec(j-1)
	end do
	endvec(nproc) = ngrid
	do j  = nproc-1,1,-1
		endvec(j) = endvec(j+1) - countvec(j+1)
	end do

	mystart = startvec(rank+1)
	myend   = endvec(rank+1)
	mycount = countvec(rank+1)

	displacement = 0
	do j = 1,nproc
		acount(j) = ncell*countvec(j)
	end do
	do j = 2,nproc
		displacement(j) = acount(j-1) + displacement(j-1)
	end do


end subroutine mpicounter


subroutine lineva1(xx,yy,x,y)


    real(8), intent(in)  :: xx(:), yy(:), x
    real(8), intent(out) :: y !, dy
    integer ix
    real(8) wx

	! 1-dim linear interpolation
	call findiw(x,xx,ix,wx)
    y  = wx*yy(ix) + (1.0d0-wx)*yy(ix+1)


end subroutine lineva1


subroutine lineva2(xgrid,ygrid,zmat,x,y,z)


    real(8), intent(in)  :: xgrid(:), ygrid(:), zmat(:,:), x, y
    real(8), intent(out) :: z
    integer ix, iy
    real(8) wx, wy, z1, z2

	! 2-dim linear interpolation
    call findiw(x,xgrid,ix,wx)
    call findiw(y,ygrid,iy,wy)

    z  = wx*wy*zmat(ix,iy) + wx*(1.0d0-wy)*zmat(ix,iy+1) + &
        (1.0d0-wx)*wy*zmat(ix+1,iy) + (1.0d0-wx)*(1.0d0-wy)*zmat(ix+1,iy+1)


end subroutine lineva2


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


subroutine tauchen(N,mu,rho,sigma,m,Z,Zprob)


    ! declare variables
    integer, intent(in) :: N
    real(8), intent(in) :: mu, rho, sigma, m
    real(8), intent(out), dimension(N) :: Z
    real(8), intent(out), dimension(N,N) :: Zprob

    real(8) c, w
    integer i, j, k

    c = (1.0d0-rho)*mu

    Z(N) = m*sqrt(sigma**2/(1.0d0-rho**2))
    Z(1) = -Z(N)
    w = (Z(N)-Z(1))/(N-1)

    do i = 2,(N-1)
        Z(i) = Z(1) + w*(i-1)
    end do

    do i = 1,N
        Z(i) = Z(i) + mu
    end do

    do j = 1,N
        do k = 1,N
            if (k==1) then
                Zprob(j,k) = cdf_normal((Z(1) - c - rho*Z(j) + w/2.0d0)/sigma)
            else if (k==N) then
                Zprob(j,k) = 1.0d0 - cdf_normal((Z(N) - c - rho*Z(j) - w/2.0d0)/sigma)
            else
                Zprob(j,k) = cdf_normal((Z(k) - c - rho*Z(j) + w/2.0d0)/sigma) - cdf_normal((Z(k) - c - rho*Z(j) - w/2.0d0)/sigma)
            end if
        end do
    end do

    return


end subroutine tauchen


real(8) function cdf_normal(x)


    real(8), intent(in) :: x


    cdf_normal = 0.5d0 * erfc(-x/sqrt(2.0d0))


end function cdf_normal


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


function gridlookup(x0,xgrid) result(ix)


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


end function gridlookup


Subroutine gridlookup2(r,grid,valplace,iloc)
! Looking for index of grid about a value
! Code from Aubhik's matlab code
! R(input) :  number of inner grids
! grid(input) : grid we have
! valplace(input) : value we have
! iloc(output) : index number
    Implicit none
    integer, intent(in) :: r
    Real(8), intent(in) :: valplace, grid(r+2)
    integer, intent(out) :: iloc
    integer :: ilow, ihigh, distance, inow
    Real(8) :: valnow

    ilow = 1
    ihigh = r+2
    distance = 2

	Do while (distance > 1)

	    inow = floor((ilow + ihigh)/2.)
	    valnow = grid(inow)

	    if (valnow > valplace) then
	        ihigh = inow
	    else
	        ilow = inow
	    end if

	    distance = ihigh - ilow

	end do

	iloc = ilow

end subroutine gridlookup2


function lini(t,ft,z) result(fz)


    real(8), intent(in) :: t(:), ft(:), z
    real(8) w, fz
    integer it, nt, ind

    nt = size(t)

    if (z<=t(1)) then
        fz = ft(1)
    elseif (z>=t(Nt)) then
        fz = ft(Nt)
    else ! t(1) < z < t(nt)

        ! ind = sum(t<z)
        ind = 0
        do it = 1,nt
            if (t(it)>=z) exit
            ind = ind+1
        end do

!        ind = sum(t<z)
        w = (t(ind+1)-z)/(t(ind+1)-t(ind))
        fz = w*ft(ind) + (1.0d0-w)*ft(ind+1)

    end if


end function lini


subroutine ginicoef(dist,grid,ngrid,ns,gini,share,negW)

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

    integer, intent(in) :: ngrid, ns
    real(8), intent(in) :: dist(:), grid(:)
    real(8), intent(out) :: gini, negW
    real(8), intent(out), dimension(ns) :: share

    integer :: i1, j1
    real(8) :: aaa, bbb, ccc, ddd
    real(8), dimension(ngrid) :: xaxis, a1, yaxis
    real(8), dimension(ngrid-1) :: dc, db
    real(8), dimension(ns) :: sharecum

    ddd = ns

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

    sharecum=0.0d0
    do j1 = 1,ns
        ccc = j1/ddd
        sharecum(j1) = lini(xaxis,yaxis,ccc)
    end do

    share = 0.0d0
    share(1) = sharecum(1)
    do j1 = 1,ns-2
        share(j1+1) = sharecum(j1+1) - sharecum(j1)
    end do
    share(ns) = 1 - sharecum(ns-1)

    ! share of nonpositive holders
    negW = 0.0d0
    do i1 = 1,ngrid
        if (grid(i1) <= 0) then
            negW = negW + dist(i1)
        end if
    end do

end subroutine ginicoef


subroutine markovss(pij,sizex,tol,xss)

    integer, intent(in) :: sizex
    real(8), intent(in), dimension(sizex,sizex) :: pij
    real(8), intent(in) :: tol
    real(8), intent(out), dimension(sizex) :: xss

    integer :: ix, jx
    real(8) :: dist
    real(8), dimension(sizex) :: xold, xnew

    xold = (1.0d0/dble(sizex))
    xnew = xold
    dist = 1000.0d0

    do while (dist > tol*2)
        do jx = 1,sizex
            xnew(jx) = 0.0d0
            do ix = 1,sizex
                xnew(jx) = xnew(jx) + pij(ix,jx)*xold(ix)
            end do
        end do
        dist = maxval(abs(log(xnew)-log(xold)))
        xold = xnew
    end do
    xss = xnew

end subroutine markovss


function std(y) result(s)


    real(8), intent(in) :: y(:)
    real(8) y1, y2, s
    integer n

    n  = size(y)
    y1 = sum(y,1)/n ! mean
    y2 = sum((y-y1)**2,1)/n ! var
    s  = y2**0.5d0


end function std


function inv(A) result(Ainv)

    real(8), dimension(:,:), intent(in) :: A
    real(8), dimension(size(A,1),size(A,2)) :: Ainv

    real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
        stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
        stop 'Matrix inversion failed!'
    end if

end function inv


end module mod_utils

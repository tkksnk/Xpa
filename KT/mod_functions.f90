module mod_functions
! my subroutines used in several projects
! January 2012, Takeki Sunakawa


implicit none


contains


function calcmu(PIe) result(mue)


    ! use mod_parameters
    real(8), intent(in) :: PIe(:,:)
    real(8), allocatable :: Pemat(:,:), mue(:)
    integer ie, ne, i

    ne = size(PIe,1)
    allocate(Pemat(ne,ne),mue(ne))

    ! stationary distribution of e
    do ie = 1,ne
        Pemat(ie,ie) = 1.0d0
    end do

    do i = 1,2000
        Pemat = matmul(Pemat,PIe)
    end do

    mue = reshape(Pemat(1,:), (/ne/))
    mue = mue/sum(mue)


end function calcmu


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

    !% symmetrically and evenly spaced between [-epsilon, epsilon] with h elements.
    !% When p = q, then then variance of shock is epsilon^2/(h-1).
    zvar = (sigmas**2)/(1.0d0 - rho**2)
    epsilon = sqrt((znum - 1.0d0)*zvar)

    Z = linspace(-epsilon, epsilon, znum)


end subroutine rouwenhorst


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


! function hpfilter(y,lam) result(yf)
! ! based on hpdata.m by Larry Christiano.  This program applies the hp filter to data
! ! Requires Intel MKL Library to solve for Ax = b with sparce matrix A
!
!
!     real(8), intent(in) :: y(:,:)
!     integer, intent(in) :: lam
!     real(8), allocatable :: a(:,:), ab(:,:), yt(:,:), yf(:,:)
!     integer i, j, d, nrhs, kl, ku
!     integer INFO
!     integer, allocatable :: IPIV(:)
!     real(8), allocatable :: WORK(:)
!
!
!     d = size(y,1)
!     nrhs = size(y,2)
!     allocate(a(d,d),yt(d,nrhs),yf(d,nrhs))
!     a = 0.0d0
!
!     do i = 3,d-2
! 	    a(i,i) = 6.0d0*lam+1.0d0
! 	    a(i,i+1) = -4.0d0*lam
! 	    a(i,i+2) = lam
! 	    a(i,i-2) = lam
! 	    a(i,i-1) = -4.0d0*lam
!     end do
!
!     a(2,2) = 1.0d0+5.0d0*lam
!     a(2,3) = -4.0d0*lam
!     a(2,4) = lam
!     a(2,1) = -2.0d0*lam
!     a(1,1) = 1.0d0+lam
!     a(1,2) = -2.0d0*lam
!     a(1,3) = lam
!
!     a(d-1,d-1) = 1.0d0+5.0d0*lam
!     a(d-1,d-2) = -4.0d0*lam
!     a(d-1,d-3) = lam
!     a(d-1,d) = -2.0d0*lam
!     a(d,d) = 1.0d0+lam
!     a(d,d-1) = -2.0d0*lam
!     a(d,d-2) = lam
!
!     ! solve a*yt = t for yt
!     allocate(IPIV(d),WORK(d))
!
!     ! standard matrix
!     !call dgetrf(d,d,a,d,IPIV,INFO)
!     !call dgetri(d,a,d,IPIV,WORK,d,INFO)
!     !yt(:,:) = matmul(a,y)
!
!     ! band martix
!     kl = 2
!     ku = 2
!     allocate(ab(2*kl+ku+1,d))
!     do j = 1,d
!         do i = max(1,j-ku),min(d,j+kl)
!             AB(kl+ku+1+i-j,j) = A(i,j)
!         end do
!     end do
!
!     yt = y
!     call dgbsv(d,kl,ku,nrhs,ab,2*kl+ku+1,IPIV,yt,d,INFO)
!
!     yf = y-yt
!
!
! end function hpfilter


function std(y) result(s)


    real(8), intent(in) :: y(:)
    real(8) y1, y2, s
    integer n

    n  = size(y)
    y1 = sum(y,1)/n ! mean
    y2 = sum((y-y1)**2,1)/n ! var
    s  = y2**0.5d0


end function std


function corr(a1,a2) result(c1)


    real(8), intent(in) :: a1(:), a2(:)
    real(8), allocatable :: b1(:), b2(:)
    real(8) c1
    integer n

    n = size(a1,1)

    allocate(b1(n),b2(n))

    b1 = a1-sum(a1,1)/n
    b2 = a2-sum(a2,1)/n
    c1 = sum(b1*b2)/(sum(b1*b1)*sum(b2*b2))**0.5d0


end function corr


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


function cov(xmat) result(S)


    real(8), intent(in) :: xmat(:,:)
    real(8), allocatable :: S(:,:), xmean(:), ymat(:,:)
    integer i, j, n, m

    n = size(xmat,1)
    m = size(xmat,2)
    allocate(S(m,m),xmean(m),ymat(n,m))

    do i = 1,m

        xmean(i) = sum(xmat(:,i),1)/n

        do j = 1,n

            ymat(j,i) = xmat(j,i) - xmean(i)

        end do

    end do

    S = matmul(transpose(ymat),ymat)/(n-1)


end function cov


function randn() result(n)
!----------------------------------------------------------------------------
! Returns a normally distributed deviate with zero mean and unit variance
! The routine uses the Box-Muller transformation of uniform deviates.
!----------------------------------------------------------------------------

    implicit none

    real(8) r, x, y, n

    do

        call random_number(x)
        call random_number(y)
        x = 2.0d0*x - 1.0d0
        y = 2.0d0*y - 1.0d0
        r = x**2 + y**2
        if (r<1.0d0) exit

    end do

    n = x*sqrt(-2.0d0*log(r)/r)


end function randn


end module mod_functions

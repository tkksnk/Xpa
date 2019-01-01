include 'mkl_vsl.fi'

module mod_randnumber


use mkl_vsl
use mkl_vsl_type

implicit none


contains


function randsample2(stream,n,w) result(y)

    ! random sampling with replacement based on randsample.m in MATLAB
    ! Written by Takeki Sunakawa
    ! August 25 2017 @ Sogang University
    !
    ! Inputs
    ! stream : Intel MKL's stream
    ! n : the length of the random variable vector
    ! w : mx1 vector of the distribution measure (pdf)
    ! Output
    ! y : nx1 vector of indices specifing the values drawn from the distribution measure

    type(vsl_stream_state), intent(in) :: stream
    integer, parameter :: method = VSL_RNG_METHOD_UNIFORM_STD
    integer, intent(in) :: n
    real(8), intent(in) :: w(:)
    real(8) sumw, cumsum
    real(8), allocatable :: x(:), p(:), edges(:)
    integer, allocatable :: y(:)
    integer m, i, j, index
    integer errcode


    ! n = size(x,1)
    m = size(w,1)
    allocate(p(m),edges(m+1),x(n),y(n))
    errcode = vdrnguniform(method,stream,n,x,0.0d0,1.0d0)

    sumw = sum(w,1)
    p = w/sumw

    edges(1) = 0.0d0
    ! cumulative sum
    cumsum = 0.0d0
    do i = 2,m+1
        cumsum = cumsum + p(i-1)
        edges(i) = cumsum
    end do

    do i=1,n

        index = 0
        do j=1,m+1

            if (edges(j)<=x(i)) then
                index = index + 1
            else
                exit
            end if

        end do

        y(i) = index

    end do

end function randsample2


function genzvec(simTT,RHOZ,SDINOVZ) result(zvec)


    integer, intent(in) :: simTT
    real(8), intent(in) :: RHOZ, SDINOVZ
    integer time
    real(8) znow
    real(8), allocatable :: zvec(:), randn(:)

    ! for mkl_vsl
    integer :: seed = 1848
    integer :: brng = vsl_brng_mt19937
    ! integer :: methodu = VSL_RNG_METHOD_UNIFORM_STD
    integer :: methodn = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    type(vsl_stream_state) :: stream
    integer :: errcode

    allocate(zvec(simTT),randn(simTT-1))

    errcode = vslnewstream(stream, brng, seed)
    ! errcode = vdrnguniform(methodu, stream, simTT-1, randn, 0.0d0, 1.0d0)
    errcode = vdrnggaussian(methodn, stream, simTT-1, randn, 0.0d0, 1.0d0)

    zvec(1) = 1.0d0
    do time = 1,simTT-1

        znow = zvec(time)
        zvec(time+1) = exp(RHOZ*log(znow) + SDINOVZ*randn(time))

    end do


end function genzvec


end module mod_randnumber

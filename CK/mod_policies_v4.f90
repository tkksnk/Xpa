module mod_policies


implicit none


contains


subroutine taxfunc(earns,capital,interest,yss,taul,a0,a1,a2,tax)

    real(8), intent(in) :: earns, capital, interest, yss, taul, a0, a1, a2
    real(8), intent(out) :: tax
    real(8) :: income, capinc, aaa, bbb, ccc ,ddd, deinc

    capinc=capital*interest
    income=earns+capinc
    deinc=income/yss
    aaa=income**(-a1)
    bbb=aaa+a2
    ccc=bbb**(-(1.0d0/a1))
    ddd=income-ccc
    !tax=a0*ddd+taul*earns    ! Total tax is the sum of progressive tax, proportioanl labor tax, and proportional capital tax
    tax=taul*earns

end subroutine taxfunc


subroutine transfers(xnow,anow,interest,wage,work,output,Trscale,Transprog,T0,trans)
    real(8), intent(in) :: xnow, anow, interest, output, Trscale, Transprog, T0, wage, work
    real(8), intent(out) :: trans
    real(8) :: income, aaaa, bbb, ccc ,ddd, capinc, deinc

    if (anow > 0) then
        income = wage*xnow*work + interest*anow
    else
        income = wage*xnow*work
    end if

    !trans=Trscale/(1.0d0+Transprog*exp(income/output))
    !trans= T0 + Trscale*((1.0d0+(income/output))**(-Transprog) )
    !trans= T0 + Trscale*((1.0d0+(income))**(-Transprog) )

    !trans= Trscale/( 1.0d0+Transprog*exp(T0*log(income)) )

	!trans= Trscale/( 1.0d0+Transprog*exp(T0*log(xnow)) )
	!trans= Trscale/( 1.0d0+exp(Transprog*log(xnow)) )
	! trans= Trscale/( 1.0d0+T0*(xnow**Transprog))
	!trans= Trscale/( 1.0d0+Transprog*exp(xnow) )

	!trans= Trscale/( 1.0d0+exp(Transprog*log(xnow)) )
    !trans= Trscale/( 1.0d0+Transprog*(xnow) )
    ! NOTE: August 2018 version
    trans= T0 + Trscale*((1.0d0+(xnow))**(-Transprog) )

    !trans= T0 + Trscale*((1.0d0+(xnow))**(-Transprog))
	!trans= T0 + Trscale*exp(-Transprog*(xnow))

end subroutine transfers


subroutine sortdist(nsize, gridold, muold, gridold2, gridnew, munew, gridnew2)
    integer, intent(in) :: nsize
    real(8), intent(in) :: gridold(nsize), muold(nsize), gridold2(nsize)
    real(8), intent(out) :: gridnew(nsize), munew(nsize), gridnew2(nsize)
    integer :: irow, krow
    real(8) :: gridtemp, mutemp, gridtemp2

    gridnew = gridold
    munew = muold
    gridnew2 = gridold2

    do irow = 1,nsize
        krow = minloc( gridnew(irow:nsize), dim=1 ) + irow - 1
        ! gridold is the sorting basis

        gridtemp = gridnew(irow)
        mutemp = munew(irow)
        gridtemp2 = gridnew2(irow)

        gridnew(irow) = gridnew(krow)
        munew(irow) = munew(krow)
        gridnew2(irow) = gridnew2(krow)

        gridnew(krow) = gridtemp
        munew(krow) = mutemp
        gridnew2(krow) = gridtemp2

    end do

end subroutine sortdist


end module mod_policies

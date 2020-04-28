SUBROUTINE labor3 (X, F, N)
!This subroutine contains the maximization problem in labor for a married couple when both spouses work
USE Model_Parameters
USE glob0
USE Utilities
INTEGER N
real(8) X(N), F(N)
real(8) :: dum2,Y1
   
    if(ind==1) then
        !Y1 is consumption to be used in the FOC
        Y1=dum
        X(1)=max(X(1),0.0001d0)
        F(1)=wagef*(1d0-(tax_labor(wagem*1d0+wagef*X(1))+tSS_employee(wagem*1d0+wagef*X(1))))-(1d0/AE)*wagef*(wagef*X(1)+wagem*1d0)*theta(1)*theta(2)*((wagem*1d0+wagef*X(1))*(1d0/AE))**(-1d0-theta(2))-chif*(X(1)**etaf)*Y1*(1d0+tc)
    else
        !Y1 is consumption to be used in the FOC
        Y1=dum
        X(1)=max(X(1),0.0001d0)
        F(1)=wagem*(1d0-(tax_labor(wagem*X(1)+wagef*1d0)+tSS_employee(wagem*X(1)+wagef*1d0)))-(1d0/AE)*wagem*(wagem*X(1)+wagef*1d0)*theta(1)*theta(2)*((wagem*X(1)+wagef*1d0)*(1d0/AE))**(-1d0-theta(2))-chim*(X(1)**etam)*Y1*(1d0+tc)
    end if
    
end subroutine labor3


subroutine labor3_hybrd(n, x, fvec, iflag)
    integer :: n, iflag
    real(8) :: x(n), fvec(n)

    call labor3(x, fvec, n)
    
end subroutine labor3_hybrd

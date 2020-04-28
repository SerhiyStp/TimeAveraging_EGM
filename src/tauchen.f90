MODULE TAUCHEN

    IMPLICIT NONE

Contains

    Subroutine tauchen_hans(STDE,RHO,ZZ,SZ,PZEE,pi)
        !=============================================================
        !DISCRETIZE A NORMAL DISTRIBUTION BY USING THE TAUCHEN(1986)
        !==============================================================
        !	Input:  STDW(standard deviation)
        !		    Log(z_t+1)=RHO*Log(z_t)+e_t+1; e~N(0,STDE^2).
        !			M is the number of st. dev. within the distribution is approximated
        !			ZZ is the number of shocks
        !	Ouput:  SZ is the set of shock values,
        !	        PZEE is the transition matrix of z
        !		    pi is a stationary distribution of z

        Implicit NONE

        Real(8), Intent(In) :: STDE,RHO
        Integer, Intent(In) :: ZZ
        !Real(8), Intent(Out) :: SZ(ZZ),PZEE(ZZ,ZZ),pi(ZZ)
        Real(8) :: SZ(ZZ),PZEE(ZZ,ZZ),pi(ZZ)
        INTEGER i,j,k
        Integer, parameter :: draws_number=1000000
        Real(8)  STDW, step, M, temp1,dum
        Real(8), Dimension(draws_number) :: normal_draws
        real(8) :: cdf1, cdf2, tmp

        M=3.0

        !==================================
        !USE PZW TO SET THE VALUE OF SHOCKS
        !==================================

        STDW=STDE/sqrt(1-RHO*RHO)
        SZ(ZZ)=M*STDW

        SZ(1)=-SZ(ZZ)
        step=2*STDW*M/FLOAT(ZZ-1)
        DO i=2,ZZ-1
            SZ(i)=SZ(i-1)+step
        Enddo

        dum=(5D-1)*step

        do i=1,ZZ
            !PZEE(i,1)=D_ANORDF((SZ(1)-RHO*SZ(i)+dum)/STDE)
            call normal_01_cdf((SZ(1)-RHO*SZ(i)+dum)/STDE, cdf2)
            PZEE(i,1) = cdf2
        end do

        do i=2,ZZ-1
            do j=1,ZZ
                !PZEE(j,i)=D_ANORDF((SZ(i)-RHO*SZ(j)+dum)/STDE)-D_ANORDF((SZ(i)-RHO*SZ(j)-dum)/STDE)
                call normal_01_cdf((SZ(i)-RHO*SZ(j)+dum)/STDE, cdf2)
                call normal_01_cdf((SZ(i)-RHO*SZ(j)-dum)/STDE, cdf1)
                PZEE(j,i) = cdf2 - cdf1
            end do
        end do

        do i=1,ZZ
            call normal_01_cdf((SZ(ZZ)-RHO*SZ(i)-dum)/STDE, cdf1)
            PZEE(i,ZZ) = 1.0d0 - cdf1
            !PZEE(i,ZZ)=1D0-D_ANORDF((SZ(ZZ)-RHO*SZ(i)-dum)/STDE)
        end do

        !=============
        !Find the stationary distribution
        !=============
        pi=1d0/ZZ
        Do i=1,100
            pi=Matmul(pi,PZEE)
        Enddo

        RETURN

    END Subroutine tauchen_hans

    subroutine normal_01_cdf ( x, cdf )

        !*****************************************************************************80
        !
        !! NORMAL_01_CDF evaluates the Normal 01 CDF.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    10 February 1999
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    AG Adams,
        !    Algorithm 39,
        !    Areas Under the Normal Curve,
        !    Computer Journal,
        !    Volume 12, pages 197-198, 1969.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) X, the argument of the CDF.
        !
        !    Output, real ( kind = 8 ) CDF, the value of the CDF.
        !
        implicit none

        real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
        real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
        real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
        real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
        real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
        real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
        real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
        real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
        real ( kind = 8 ), parameter :: b1 = 3.8052D-08
        real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
        real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
        real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
        real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
        real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
        real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
        real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
        real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
        real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
        real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
        real ( kind = 8 ) cdf
        real ( kind = 8 ) q
        real ( kind = 8 ) x
        real ( kind = 8 ) y
        !
        !  |X| <= 1.28.
        !
        if ( abs ( x ) <= 1.28D+00 ) then

            y = 0.5D+00 * x * x

            q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
                + a6 / ( y + a7 ) ) ) )
            !
            !  1.28 < |X| <= 12.7
            !
        else if ( abs ( x ) <= 12.7D+00 ) then

            y = 0.5D+00 * x * x

            q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
                + b2 / ( abs ( x ) + b3 &
                + b4 / ( abs ( x ) - b5 &
                + b6 / ( abs ( x ) + b7 &
                - b8 / ( abs ( x ) + b9 &
                + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
            !
            !  12.7 < |X|
            !
        else

            q = 0.0D+00

        end if
        !
        !  Take account of negative X.
        !
        if ( x < 0.0D+00 ) then
            cdf = q
        else
            cdf = 1.0D+00 - q
        end if

        return
    end subroutine normal_01_cdf

END MODULE TAUCHEN

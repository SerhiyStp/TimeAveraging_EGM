module MyInterpolation
    implicit none
    
contains

    function S_QDVAL(xval, xdata, ydata, ndata) result(yval)
        real(8) :: xval
        integer :: ndata
        real(8) :: xdata(ndata), ydata(ndata)
        integer :: left
        real(8) :: yval
        real(8) :: t1, t2, t3
        real(8) :: y1, y2, y3
        real(8) :: dif1, dif2
        
        left = 1
        call r8vec_bracket3 ( ndata, xdata, xval, left )    
        left = min(ndata - 2, left)
        
        t1 = xdata(left)
        t2 = xdata(left+1)
        t3 = xdata(left+2)        
        
        y1 = ydata(left)
        y2 = ydata(left+1)
        y3 = ydata(left+2)

        dif1 = ( y2 - y1 ) / ( t2 - t1 )
        dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
             - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

        yval = y1 + ( xval - t1 ) * ( dif1 + ( xval - t2 ) * dif2 )        
    
    end function S_QDVAL
    
end module MyInterpolation

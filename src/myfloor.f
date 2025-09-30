      INTEGER FUNCTION myfloor( x )
* This function rounds  x  towards -infinity

      IMPLICIT NONE
      DOUBLE PRECISION  x

      if ( x .GT. 0d00 ) then
         myfloor = int( x )
      else
         myfloor = int( x ) - 1
      endif

      RETURN
      END

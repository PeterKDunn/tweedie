      INTEGER FUNCTION ceil( x )
* This function rounds  x  towards +infinity

      IMPLICIT NONE
      DOUBLE PRECISION  x

      if ( x .GT. 0d00 ) then
         ceil = int( x ) + 1
      else
         ceil = int( x )
      endif

      RETURN
      END
      
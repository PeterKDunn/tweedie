      SUBROUTINE findZeroSmallp(t, f, df)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, t
      DOUBLE PRECISION f, df, Imk, Imdk
      INTEGER m
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      COMMON /mparam/ m 

      CALL findImk(t, Imk) 
      CALL findImkd(t, Imdk)
      
      f  = sin(Imk + (t * Cy) ) - sin( Imk)
      df = (Imdk + Cy) * cos(Imk + (t * Cy) ) -
     &       Imdk * cos(Imk)
*      write(*,*) "f = ", f
*      write(*,*) "df = ", df
      RETURN 
      END
      
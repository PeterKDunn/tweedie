      SUBROUTINE findImkM(t, f, df)

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, t, pi, f, df
      DOUBLE PRECISION Imk, omega, pindex, front, alpha
      DOUBLE PRECISION Imdk
      INTEGER m
      COMMON /params/ Cp, Cy, Cmu, Cphi
      COMMON /mparam/ m 

      pi = 4.0d0 * DATAN(1.0d0)

      CALL findImk(t, Imk) 
      CALL findImkd(t, Imdk)
      
      f  = Imk - DBLE(m) * Cy
      df = Imdk

      RETURN 
      END
      
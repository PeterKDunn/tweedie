      SUBROUTINE findImdkZero(t, f, df)
*     Evaluates Im(k'(t)) and Im(k''(t)), for finding kmax
*     Returns both.

      IMPLICIT NONE
      DOUBLE PRECISION t, f, df, Imdk, Imddk
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      CALL findImkd( t, Imdk ) 
      CALL findImkdd(t, Imddk)
*      write(*,*) '+findImdk', Cp, Cy, Cmu, Cphi, t
*      write(*,*) '+findImdk', Imdk, Imddk

      
      f  = Imdk
      df = Imddk

      RETURN
      END

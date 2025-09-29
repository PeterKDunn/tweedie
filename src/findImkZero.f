      SUBROUTINE findImkZero(t, f, df)
*     Evaluates Im(k'(t)) and Im(k''(t)), for finding kmax
*     Returns both.

      DOUBLE PRECISION t, f, df, Imdk, Imddk
      DOUBLE PRECISION Cp, Cy, Cmu, Cphi
      COMMON /params/ Cp, Cy, Cmu, Cphi

      CALL findImkd(t, Imdk) 
      CALL findImkdd(t, Imddk)
      
      f  = Imk
      df = Imdk

      RETURN
      END

      DOUBLE PRECISION FUNCTION findKmaxSP()
      
*     Compute a starting guess for t in
*       mu * cos(atan((1-p)*t*phi/mu^(1-p)) / (1-p)) / 
*       cos(atan((1-p)*t*phi/mu^(1-p)))**(1/(1-p)) - y = 0
*
*     Works for p > 1; smoothly interpolates small/large t

      DOUBLE PRECISION Cmu, Cy, Cphi, Cp, pi
      DOUBLE PRECISION tsmall, tlarge, abs1mp
      DOUBLE PRECISION omegaInf, slope
      DOUBLE PRECISION slope1, slopesmall, slopelarge
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      abs1mp = DABS(1.0d00 - Cp)
      pi = 4.0d00 * DATAN(1.0d00)

*     We find a small-t approx, a large-t approx, and a combined approx.
*     The SP should be the FIRST of these whose slope is NEGATIVE

*     Small-t approximation
* SEEMS TO WORK OK for p > 2 (and for all y > mu).
* BUT NOT FOR 1<p<2, and y < mu
*      tsmall = DSQRT(2.0d00 * (Cmu - Cy)/Cmu) * 
*     &         Cmu**(1.0d00 - Cp) / Cphi
 
       write(*,*) "Update approach in text: largest where slope -ive"
 
       omegaInf = (pi / 2.0d00) *  
     &            (1.0d00 - Cp)/(2.0d00*Cp - 1.0d00)
       tsmall = Cmu**(1.0d00 - Cp) / ( (1.0d00 - Cp)) * 
     &          DTAN(omegaInf)
       CALL findImkd(tsmall, slope)

      IF (slope .LE. 0.0d0) THEN
        findKmaxSP = tsmall

        RETURN
      ENDIF
      
*      Large-t approximation
       tlarge = (Cmu / Cy)**(Cp - 1.0d00) * 
     &          Cmu**(1.0d00 - Cp) / (Cphi * abs1mp)

      CALL findImkd(tlarge, slope)

      IF (slope .LE. 0.0d0) THEN
        findKmaxSP = tsmall

        RETURN
      ENDIF

*     Smooth interpolation: sum of small + large contributions
      findKmaxSP = tsmall + tlarge
      write(*,*) "FIX: If slope at tlarge not negative, need to 
     &            make t larger, not compromise"
     
      RETURN
      END

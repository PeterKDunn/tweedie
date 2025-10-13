      DOUBLE PRECISION FUNCTION findKmaxSP()
      
*     Compute a starting guess for t in
*       mu * cos(atan((1-p)*t*phi/mu^(1-p)) / (1-p)) / 
*       cos(atan((1-p)*t*phi/mu^(1-p)))**(1/(1-p)) - y = 0
*
*     Works for p > 1; smoothly interpolates small/large t

      IMPLICIT NONE
      DOUBLE PRECISION Cmu, Cy, Cphi, Cp, pi
      DOUBLE PRECISION tsmall, tlarge, abs1mp, dslope
      DOUBLE PRECISION omegaInf, slope, candidate, front
      DOUBLE PRECISION oldCandidate
      INTEGER n
      LOGICAL pSmall, found
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      abs1mp = DABS(1.0d00 - Cp)
      pi = 4.0d00 * DATAN(1.0d00)

*     We find a small-t approx, a large-t approx, and a combined approx.
*     The SP should be the FIRST of these whose slope is NEGATIVE

       write(*,*) "Update approach in text: largest where slope -ive"
      IF (pSmall .AND. (Cy .LT. Cmu) ) THEN
*       Find the turning points.
*       Find the slope at each of these until one is negative.
*       This, it is between them...?

*     Try small-t approximation
* SEEMS TO WORK OK for p > 2 (and for all y > mu).
* BUT NOT FOR 1<p<2, and y < mu???

        tsmall = DSQRT(2.0d00 * (Cmu - Cy)/Cmu) * 
     &         Cmu**(1.0d00 - Cp) / Cphi
         findKmaxSP = tsmall
         write(*,*) ">>> tsmall", tsmall

      ELSE
        omegaInf = (pi / 2.0d00) *  
     &             (1.0d00 - Cp)/(2.0d00*Cp - 1.0d00)
        tsmall = Cmu**(1.0d00 - Cp) / ( (1.0d00 - Cp)) * 
     &           DTAN(omegaInf)
        CALL findImkd(tsmall, slope)
        RETURN
        IF (slope .LE. 0.0d0) THEN
          findKmaxSP = tsmall
      write(*,*) " findKmaxSP FOUND: tsmall, slope", tsmall, slope
          RETURN
        ENDIF
      ENDIF
*      Large-t approximation
       tlarge = (Cmu / Cy)**(Cp - 1.0d00) * 
     &          Cmu**(1.0d00 - Cp) / (Cphi * abs1mp)

      CALL findImkd(tlarge, slope)

      IF (slope .LE. 0.0d0) THEN
        findKmaxSP = tlarge
      write(*,*) "findKmaxSPFOUND: tlarge, slope", tlarge, slope

        RETURN
      ENDIF

*     Smooth interpolation: sum of small + large contributions
      findKmaxSP = tsmall + tlarge
      write(*,*) "FIX: If slope at tlarge not negative, need to 
     &            make t larger, not compromise"
     
      RETURN
      END

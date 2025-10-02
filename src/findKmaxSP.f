      DOUBLE PRECISION FUNCTION findKmaxSP()
      
*     Compute a starting guess for t in
*       mu * cos(atan((1-p)*t*phi/mu^(1-p)) / (1-p)) / 
*       cos(atan((1-p)*t*phi/mu^(1-p)))**(1/(1-p)) - y = 0
*
*     Works for p > 1; smoothly interpolates small/large t

      DOUBLE PRECISION Cmu, Cy, Cphi, Cp, pi
      DOUBLE PRECISION tsmall, tlarge, abs1mp
      DOUBLE PRECISION omegaInf
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      abs1mp = DABS(1.0d00 - Cp)
      pi = 4.0d00 * DATAN(1.0d00)

*     Small-t approximation
*      tsmall = DSQRT(2.0d00 * (Cmu - Cy)/Cmu) * 
*     &         Cmu**(1.0d00 - Cp) / Cphi
 
      omegaInf = (pi / 2.0d00) * 
     &           (1.0d00 - Cp)/(2.0d00*Cp - 1.0d00)
      tsmall = Cmu**(1.0d00 - Cp) / ( (1.0d00 - Cp)) * 
     &         DTAN(omegaInf)
      
*     Large-t approximation
      tlarge = (Cmu / Cy)**(Cp - 1.0d00) * 
     &          Cmu**(1.0d00 - Cp) / (Cphi * abs1mp)

*     Smooth interpolation: sum of small + large contributions
      findKmaxSP = tsmall + tlarge
      write(*,*) "The three starting points"
      write(*,*) tsmall, tlarge, findKmaxSP
      
      RETURN
      END

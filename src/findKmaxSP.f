      DOUBLE PRECISION FUNCTION findKmaxSP()
      
*     Compute a starting guess for t in
*       mu * cos(atan((1-p)*t*phi/mu^(1-p)) / (1-p)) / 
*       cos(atan((1-p)*t*phi/mu^(1-p)))**(1/(1-p)) - y = 0
*
*     Works for p > 1; smoothly interpolates small/large t

      DOUBLE PRECISION Cmu, Cy, Cphi, Cp
      DOUBLE PRECISION tsmall, tlarge, abs1mp
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

      abs1mp = DABS(1.0d00 - Cp)

*     Small-t approximation
      tsmall = DSQRT(2.0d00 * (Cmu - Cy)/Cmu) * 
     &         Cmu**(1.0d00 - Cp) / Cphi

*     Large-t approximation
      tlarge = (Cmu / Cy)**(Cp - 1.0d00) * 
     &          Cmu**(1.0d00 - Cp) / (Cphi * abs1mp)

*     Smooth interpolation: sum of small + large contributions
      findKmaxSP = tsmall + tlarge

      RETURN
      END

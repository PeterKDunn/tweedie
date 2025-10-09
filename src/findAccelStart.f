
      SUBROUTINE findAccelStart(tRightMost)

* For the case 1<p,2, find when to start the Sidi acceleration.
* It should start after the final turning point of Re k(t).

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, pi
      DOUBLE PRECISION omegaRM, tRightMost
      DOUBLE PRECISION lRightMostD
      INTEGER lRightMost, myfloor
      LOGICAL pSmall
      EXTERNAL myfloor
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      
      
      pi = 4.0d00 * DATAN(1.0d0)

      IF ( (Cp .GT. 1.4d00) .AND.
     &     (Cp .LT. 1.6d00) ) THEN
     
        tRightMost = 1.0d0
      ELSE
        lRightMostD = DACOS( (Cy / Cmu)**(Cp - 1.0d0) ) / 
     &               (pi * (Cp - 1.0d0) ) 
        lRightMost = myfloor(lRightMostD) + 1

*     Add one, since without it we have the TP... but it can sometimes
*     take some time to settle down

        omegaRM = lRightMost * pi * (1.0d0 - Cp)
        tRightMost = tan(omegaRM) * Cmu ** (1.0d0 - Cp) / 
     &                            ( Cphi * (1.0d0 - Cp) )
      ENDIF
      
      
      RETURN
      END
      

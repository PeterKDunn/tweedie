
      SUBROUTINE findAccelStart(tRightMost)

* For the case 1<p,2, find when to start the Sidi acceleration.
* It should start after the final turning point of Re k(t).

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, pi
      DOUBLE PRECISION lRightMost, omegaRM, tRightMost
      DOUBLE PRECISION myfloor
      LOGICAL pSmall
      EXTERNAL myfloor
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      
      
      pi = 4.0d00 * DATAN(1.0d0)

      lRightMost = DACOS( (Cy / Cmu)**(Cp - 1.0d0) ) / 
     &             (pi * (Cp - 1.0d0) ) 
      lRightMost = myfloor(lRightMost) + 1
*     Add one, since without it we have the TP... but it can sometimes
*     take some time to settle down

      omegaRM = lRightMost * pi * (1.0d0 - Cp)
      tRightMost = tan(omegaRM) * Cmu ** (1.0d0 - Cp) / 
     &                            ( Cphi * (1.0d0 - Cp) )
      write(*,*) "findAccelStart"
      write(*,*) "l, omega, t"
      write(*,*) lRightMost, omegaRM, tRightMost
     

      RETURN
      END
      

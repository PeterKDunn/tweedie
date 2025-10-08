
      SUBROUTINE findAccelStart(tRightMost)

* For the case 1<p,2, find when to start the Sidi acceleration.
* It should start after the final turning point of Re k(t).

      IMPLICIT NONE
      DOUBLE PRECISION Cp, Cmu, Cphi, Cy, pi
      DOUBLE PRECISION omegaRM, tRightMost
      DOUBLE PRECISION lRightMostD, a, r, x
      INTEGER lRightMost, myfloor
      LOGICAL pSmall
      EXTERNAL myfloor
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall
      
      
      pi = 4.0d00 * DATAN(1.0d0)

      write(*,*) "findAccelStart"
      
      IF ( (Cp .GT. 1.4d00) .AND.
     &     (Cp .LT. 1.6d00) ) THEN
     
*        a = Cphi * DSQRT(Cmu) / 2.0d0
*        r = Cy / Cmu
*        x = (-(2.0d0 * r + 1.0d0) + DSQRT(1.0d0 + 8.0d0 * r) ) /
*     &      (2.0d0 * r)
*        tRightMost = DSQRT(x) / a
*        write(*,*) "tRightmost when p=1.5:", tRightMost

        tRightMost = 1.0d0
        write(*,*) "FIX IN findAccelStart!!!"
      ELSE
        lRightMostD = DACOS( (Cy / Cmu)**(Cp - 1.0d0) ) / 
     &               (pi * (Cp - 1.0d0) ) 
      write(*,*) "l (before floor() ): ", lRightMostD
        lRightMost = myfloor(lRightMostD) + 1
      write(*,*) "l (after floor() ): ", lRightMost

*     Add one, since without it we have the TP... but it can sometimes
*     take some time to settle down

        omegaRM = lRightMost * pi * (1.0d0 - Cp)
        tRightMost = tan(omegaRM) * Cmu ** (1.0d0 - Cp) / 
     &                            ( Cphi * (1.0d0 - Cp) )
        write(*,*) "l, omega, t"
        write(*,*) lRightMost, omegaRM, tRightMost
     
      ENDIF
      
      
      RETURN
      END
      

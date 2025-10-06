      DOUBLE PRECISION FUNCTION DFintegrand(t)

*     A function to be numerically integrated in the DF.
*     Note that
*        lim (n->infty) exp( Re k) = exp( -lambda).
*     DFfun2  is used for the conditional CDF for 1 < p < 2; DFfun  for p > 2.
*
*     IN:  p, phi, y, mu, t
*     OUT: DFfun

      DOUBLE PRECISION t, Cp, Cy, Cmu, Cphi
      DOUBLE PRECISION Imk, Rek, lambda
      LOGICAL pSmall
      COMMON /params/ Cp, Cy, Cmu, Cphi, pSmall

* MAJOR VARIABLES:
*   t          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*   lambda     : P(Y = 0) = exp( -lambda ) when 1 < p < 2
 

*     Check for when t = 0, which should never actually happen 
      IF (DABS(t) .LT. 1.0d-14) THEN
*        IF ( t .EQ. 0.0d00 ) THEN
         DFintegrand = Cmu - Cy
         write(*,*) "!!!!! DFint: should never happen: t = 0 !!!!!"
      ELSE
        CALL findImk(t, Imk)
        CALL findRek(t, Rek)
        
        DFintegrand = DEXP( Rek ) * DSIN( Imk ) / t

        IF (pSmall) THEN
          CALL findLambda(lambda)
          DFintegrand = DFintegrand - 
     &                    ( DEXP(Rek) * DSIN( Imk + (t * Cy) ) ) / t
        ENDIF

      ENDIF
          write(*,*) " "
*         write(*,*) "t", t
*         write(*,*) "DFintegrand", DFintegrand
*         write(*,*) "Imk", Imk
*         write(*,*) "DSIN(Imk)", DSIN(Imk)
*         write(*,*) "Rek/t", Rek/t

      RETURN
      END
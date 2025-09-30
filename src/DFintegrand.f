      DOUBLE PRECISION FUNCTION DFintegrand(t)

*     A function to be numerically integrated in the DF.
*     Note that
*        lim (n->infty) exp( Re k) = exp( -lambda).
*     DFfun2  is used for the conditional CDF for 1 < p < 2; DFfun  for p > 2.
*
*     IN:  p, phi, y, mu, t
*     OUT: DFfun

      DOUBLE PRECISION t, Cp, Cy, Cmu, Cphi, aimrerr
      DOUBLE PRECISION calclambda, Imk, Rek, lambda
      COMMON /params/ Cp, Cy, Cmu, Cphi, aimrerr

* MAJOR VARIABLES:
*   p          : the index in the variance function, V(mu) = phi * mu^p
*   phi        : the dispersion parameter
*   y          : the point at which the DF is to be evaluated
*   mu         : the mean value
*   t          : the internal variable over which to integrate in this function
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf
*   lambda     : P(Y = 0) = exp( -lambda ) when 1 < p < 2
 


      lambda = calclambda( Cp, Cphi, Cmu )
      IF ( t .EQ. 0.0d00 ) THEN
         DFintegrand = Cmu - Cy
* ?????
      ELSE
         CALL findImk(t, Imk)
         CALL findRek(t, Rek)

         DFintegrand = DEXP( Rek ) * DSIN( Imk ) / t
*         write(*,*) " "
*         write(*,*) "t", t
*         write(*,*) "Imk", Imk
*         write(*,*) "DSIN(Imk)", DSIN(Imk)
*         write(*,*) "Rek/t", Rek/t
      ENDIF

      RETURN
      END
***
*     Author:         Peter Dunn
*     Creation date:  15 February 1999
*     Last revision:  25 Septembeer 2025
*
* IN R, print statements for debugging as follows:
*     call dblepr("The value of y is ",-1, y, 1)
*     call intpr("The value of exact is ",-1, exact, 1)
******************************************************************
** Routines common to density and DF
******************************************************************


*******************************************************************
      double precision function calclambda(p, phi, mu )
*     Calculates lambda, such that  P(Y = 0) = exp( -lambda )

      double precision  p,  phi, mu, twomp

      twomp = 2.0d00 - p      
      calclambda = (mu ** twomp) / (phi * twomp)

      return
      end

*******************************************************************
****************************************************************
      subroutine kt(p, phi, y, mu, x, rl, im)
*     Calculates real and imaginary parts of  k(t)  for the DF: k(t) = K(it) - ity
*
*     IN:  p, phi, y, mu, x
*     OUT: rl, im
*
      double precision  p, phi, y, mu, x, rl, im,
     &                  omega, alpha, front, denom

* MAJOR VARIABLES:
*   rl         : the real part of the cgf
*   im         : the imaginary part of the cgf

      omega = datan( ( 1.0d00 - p ) * x * phi /
     &               ( mu ** (1.0d00 - p) ) )

      alpha = ( 2.0d00 - p ) / ( 1.0d00 - p )
      front = ( mu ** ( 2.0d00 - p ) ) / 
     &         ( phi * ( 2.0d00 - p ) )
      denom = dcos( omega ) ** alpha

      rl = front * ( dcos ( omega * alpha ) / denom - 1)
      im = ( front * dsin ( omega * alpha ) / denom ) - 
     &         x * y

      return
      end

****************************************************************
****************************************************************
      subroutine kdasht(p, phi, y, mu, x, rl, im)
*     Calculates real and imaginary parts of  k'(t)  for the DF: k(t) = K(it) - ity
*
*     IN:  p, phi, y, mu, x
*     OUT: rl, im

      double precision  p, phi, y, mu, x, rl, im,
     &                  omega, alpha, denom

* MAJOR VARIABLES:
*   rl         : the real part of the derivative of the cgf
*   im         : the imaginary part of the derivative of the cgf

      omega = datan( ( 1.0d00 - p ) * x * phi /
     &               ( mu ** (1.0d00 - p) ) )

      alpha = 1.0d00 / ( 1.0d00 - p )
      denom = dcos( omega ) ** alpha

      rl = -mu * ( dsin( omega * alpha ) / denom )
      im =  mu *   dcos( omega * alpha ) / denom - y

      return
      end

****************************************************************
*****************************************************************
      subroutine kddasht( p, phi, mu, x, ddrl, ddim)
*     Calculates real and imaginary parts of  k''(t)  for the DF: k(t) = K(it) - ity

*     IN:   p, phi, mu, x
*     OUT:  ddrl, ddim

      double precision  p, phi, mu, x, ddrl, ddim, 
     &                  denom, alpha, omega

* MAJOR VARIABLES:
*     ddrl      : the real part of the second derivative
*     ddim      : the imaginary part of the second derivative

      omega = datan( ( 1.0d00 - p ) * x * phi /
     &               ( mu ** (1.0d00 - p) ) )

      alpha = p / ( 1.0d00 - p )

      denom = dcos( omega ) ** alpha

      ddrl = -phi * mu**p *    dsin( omega * alpha ) 
     &            / denom
      ddim = -phi * mu**p  * ( dcos( omega * alpha ) 
     &            / denom )

      return
      end

******************************************************************

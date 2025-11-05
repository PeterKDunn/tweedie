tweedie_Integrand <- function(y, power, mu, phi, t, type = "PDF", yLimits = NULL){
  
  # BEGIN: Define function to be used
  k <- function(p, mu, phi, y, t){
    front <- ( mu^(2 - p) ) / (phi * (2 - p) )
    alpha <- (2 - p)/(1 - p) 
    omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )

    Real <- front * ( cos(omega * alpha) / ( cos(omega)^alpha) - 1) 
    Imag <- front *   sin(omega * alpha) / ( cos(omega)^alpha ) - (t * y)

    return( list(Real = Real,
                 Imag = Imag) )
  }
  # END: Define function to be used

  
  kvals <- k(p = power, 
             mu = mu, 
             phi = phi,
             y = y, 
             t = t)
  k_Imag <- kvals$Imag
  k_Real <- kvals$Real
  

  
  par( mfrow = c(2,2))
  
  ### PLOT 1: Im k(t) vs t
  plot(k_Imag ~ t,
       main = expression( bold(Imaginary~part~of~italic(k)*(italic(t)))),
       xlab = expression(Values~of~italic(t)),
       ylab = expression(Im*"("*italic(k)*")"),
       las = 1,
       lwd = 2,
       type = "l")
  
  # Determine  m  values to display
  mMax <- max(0, max(kvals$Imag)) / pi + ifelse(type == "CDF", 0, pi/2)
  mMin <- min(0, min(kvals$Imag)) / pi - ifelse(type == "CDF", 0, pi/2)
  mValues <- ( floor(mMin) : floor(mMax) )

  # Adornments
  mtext(text = expression(atop(italic(m), " ") ),
        side = 4,
        las = 1,
        cex = 0.9,
        line = 2,
        adj = 1,
        at = max(k_Imag) )
  abline(h = 0, 
         col="grey")
  abline(h = mValues * pi + ifelse(type=="PDF", pi/2, 0),
         lty = 2,
         col="grey")
  axis(side = 4, 
       at = mValues * pi + ifelse(type=="PDF", pi/2, 0),
       las = 1,
       cex = 0.8,
       label = mValues )

  
  
  ### PLOT 2: sin( Im k(t) ) vs t
  plot(sin(k_Imag) ~ t,
       main = expression(sin(Im*"("*italic(k)*")")),
       xlab = expression(Values~of~italic(t)),
       ylab = expression(sin(Im*"("*italic(k)*")")),
       las = 1,
       lwd = 2,
       type = "l")
  abline(h = 0, 
         col="grey")

  
  
  
  ### PLOT 3: sin( Re k(t) ) vs t
  plot(k_Real ~ t,
       main = "Real part of k(t)",
       xlab = expression(Values~of~italic(t)),
       ylab = expression(Re*"("*italic(k)*")"),
       las = 1,
       lwd = 2,
       type = "l")
  abline(h = 0, 
         col="grey")
  
  
  ### PLOT 4: Integrand
  
  if(type == "PDF" ){
    igrand <- exp(k_Real) * cos(k_Imag)
    envelope <- exp(k_Real)
  } else {
    igrand <- exp(k_Real) * sin(k_Imag)/t
    envelope <- exp(k_Real) / t
  }
  
  plot(  igrand ~ t,
        main = "Integrand",
        xlab = expression(Values~of~italic(t)),
        ylab = "Integrand",
        las = 1,
        ylim = yLimits,
        #ylim = c(-0.0001, 0.0001),
        lwd = 2,
        type = "l")
  abline(h = 0, 
         col="grey")
  lines(x = t,
        y = envelope,
        lty = 2,
        col = "grey")
  lines(x = t,
        y = -envelope,
        lty = 2,
        col = "grey")
  
  return(list(Real = k_Real,
              Imag = k_Imag))
  
}
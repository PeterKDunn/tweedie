source("tweedieFunctions.R")

########################################################################
## TMP
y <- 0.03
mu <- 2
phi <- 1
p <- 3.5
# y <- 0.01
t <- seq(700, 2500,
         length = 1000)
########################################################################


kdashdash <- function(p, mu, phi, y, t){
  omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )
  cat("omega=", omega,"\n")
  cat("pindex=", p/(1-p),"\n")
  print( head( cbind(
                     t, 
                     omega, 
                     omegaInner = ( (1 - p) * t * phi) / mu^(1 - p),
                     front = -phi * mu ^ (p/(1 - p)),
                     TOP = sin(omega * p / (1 - p) ),
                     BOTTOM = (cos(omega)^(p/(1 - p))) ,
                     cosOmega = cos(omega),
                     df = -phi * mu ^ (p/(1 - p)) * sin(omega * p / (1 - p) ) / (cos(omega)^(p/(1 - p))) 
        ) ) )
  -phi * mu ^ (p/(1 - p)) * sin(omega * p / (1 - p) ) / (cos(omega)^(p/(1 - p))) 
  
}
kmax_est <- function(p, mu, phi, y){
  front <- (mu^(1-p))/phi
  front * (
        sqrt(2*(mu - y)/mu) + 1/(p - 1)*(mu/y)
  )
}
  



lambda <- mu^(2-p) / (phi * (2-p) )
p0 = exp(-lambda)

kd1 <- kdash(p = p, 
             mu = mu, 
             phi = phi, 
             y = y, 
             t = t)

kd2 <- kdashdash(p = p, 
                 mu = mu, 
                 phi = phi,
                 y = y, 
                 t = t)
kvals <- k(p = p, 
           mu = mu, 
           phi = phi,
           y = y, 
           t = t)

print( tail(cbind(y, 
                  t, 
                  Imag=kvals$Imag, 
                  Real=kvals$Real)))

zs <- c(10.2, 214.8, 404.7, 540.2, 659.9, 771.1, 876.8, 978.4, 1077.1)

############################################################################3

par( mfrow=c(3,2))
plot(kvals$Imag ~ t,
     main = expression( bold(Imaginary~part~of~italic(k)*(italic(t)))),
     xlab = expression(Values~of~italic(t)),
     ylab = "Im k",
     las = 1,
     lwd = 2,
     type = "l")
abline(h = 0, 
       col="grey")
abline(v = 0.03, 
       col="grey")
abline(h = (-14:14)*pi, 
       col="grey")
axis(side = 4, 
     at = c(-14:14)*pi,
     las = 1,
     cex = 0.8,
     label = (-14:14) )
abline(v = zs,
       col = "grey",
       lty = 2)



plot(sin(kvals$Imag) ~ t,
     main = "sin(Im k)",
     xlab = expression(Values~of~italic(t)),
     ylab = "sin(Im k)",
     las = 1,
     lwd = 2,
     type = "l")
abline(h = 0, 
       col="grey")
abline(v = 0.03, 
       col="grey")
abline(v = zs,
       col = "grey",
       lty = 2)





plot(kvals$Real ~ t,
     main = "Real part of k(t)",
     xlab = expression(Values~of~italic(t)),
     ylab = "Re k(t)",
     las = 1,
     #ylim = c(-0.0001, 0.0001),
     lwd = 2,
     type = "l")
abline(h = 0, 
       col="grey")



plot( exp(kvals$Real) * sin(kvals$Imag)/t ~ t,
     main = "Integrand",
     xlab = expression(Values~of~italic(t)),
     ylab = "Integrand",
     las = 1,
     #ylim = c(-0.0001, 0.0001),
     lwd = 2,
     type = "l")
abline(h = 0, 
       col="grey")
lines(x = t,
      y = exp(kvals$Real)/t,
      lty = 2,
      col = "grey")
lines(x = t,
      y = -exp(kvals$Real)/t,
      lty = 2,
      col = "grey")



plot( kd1~ t,
     main = "Derivative k'(t)",
     xlab = expression(Values~of~italic(t)),
     ylab = "deriv",
     las = 1,
     #ylim = c(-0.0001, 0.0001),
     lwd = 2,
     type = "l")
abline(h = 0, 
       col="grey")
lines(x = t,
      y = exp(kvals$Real)/t,
      lty = 2,
      col = "grey")



plot(kd2 ~ t,
     main = "Second deriv k''(t)",
     xlab = expression(Values~of~italic(t)),
     ylab = "2nd deriv",
     las = 1,
     #ylim = c(-0.0001, 0.0001),
     lwd = 2,
     type = "l")
abline(h = 0, 
       col="grey")
lines(x = t,
      y = exp(kvals$Real)/t,
      lty = 2,
      col = "grey")




if ( p < 2 ){
  cat("lambda = ", lambda, " and P(Y=0) = exp(-lambda) = ", exp(-lambda), "\n") 
}
if (mu > y ){
  cat("ESTIMATE of k_max:", kmax_est(p, mu, phi, y), "\n" )
}

# library(tweedie); ptweedie.inversion(0.05, mu=1.4, phi=0.74, power=3, exact=TRUE)
# y /q = 0.1 (or was it 0.01>???) DIES. 
# kmax must be huge, as zeros are sought for m = 390 000 and greater...
# but...that seems wrong


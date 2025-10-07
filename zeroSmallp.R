
## TMP
mu <- 1
phi <- 4
p <- 1.5
# y <- 0.01
y <- 3
t <- seq(0, 3,
         length = 1000)



kvals <- k(p = p, 
                mu = mu, 
                phi = phi,
                y = y, 
                t = t)


zeroFn <- function(p, mu, phi, y, t){
  rk <- k(p = p,
          mu = mu, 
          phi = phi, 
          y = y, 
          t = t)
  sin( rk$Imag + t * y ) - sin( rk$Imag)
}

plot( zeroFn(p, mu, phi, y, t) ~ t,
      type = "l")
abline(h = 0,
       col = "grey")



n = (1:20)
abline(v = n * pi / y,
       col = "grey",
       lty = 2)


###

rk <- k(p = p,
        mu = mu, 
        phi = phi, 
        y = y, 
        t = t)



plot( (rk$Real * (sin(rk$Imag + t*y) - sin(rk$Imag) ) /t) ~ t,
      type = "l")
abline(h = 0,
       col = "grey")

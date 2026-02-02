## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tweedie)

## ----TWexample----------------------------------------------------------------
library(statmod)

set.seed(96)
N <- 25
discrete_events <- rpois(N, lambda = 1)
y <- vector("numeric", N)

for (i in 1:N){
  y[i] <- ifelse(discrete_events[i] == 0,
                 0,
                 sum( rgamma(discrete_events, shape = 1, rate = 1)) )
}

mod.tw <- glm(y ~ 1, 
              family = statmod::tweedie(var.power = 1.5, link.power = 0) )
              # link.power = 0  means the log-link

## ----TWrandom-----------------------------------------------------------------
tweedie::rtweedie(10, xi = 1.1, mu = 2, phi = 1)

## ----TWplots------------------------------------------------------------------
y <- seq(0, 2, length = 100)
xi <- 1.1
mu <- 0.5
phi <- 0.4

twden <- tweedie::dtweedie(y, xi = xi, mu = mu, phi = phi)  
twdtn <- tweedie::ptweedie(y, xi = xi, mu = mu, phi = phi)

par(mfrow = c(1, 2))
plot( twden[y > 0] ~ y[y > 0], 
      type ="l",
      lwd = 2,
      xlab = expression(italic(y)),
      ylab = "Density function")
points(twden[y==0] ~ y[y == 0],
      lwd = 2,
      pch = 19,
      xlab = expression(italic(y)),
      ylab = "Distribution function")

plot(twdtn ~ y,
     type = "l",
      lwd = 2,
     ylim = c(0, 1),
      xlab = expression(italic(y)),
      ylab = "Distribution function")

## ----TWplots2-----------------------------------------------------------------
par(mfrow = c(1, 2))
tweedie::tweedie_plot(y, xi = xi, mu = mu, phi = phi,
                      ylab = "Density function",
                      lwd = 2)
tweedie::tweedie_plot(y, xi = xi, mu = mu, phi = phi, 
                      ylab = "Distribution function",
                      lwd = 2, 
                      ylim = c(0, 1), 
                      type = "cdf")

## ----TWqqplot-----------------------------------------------------------------
library(tweedie)

qqnorm( statmod::qresid(mod.tw) )


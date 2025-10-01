# library(tweedie); ptweedie.inversion(1, mu=1.4, phi=0.74, power=3, exact=TRUE)
source("./tweedieFunctions.R")

y=1 
mu=1.4
phi=0.74
p=3

## INITIAL
area0 = integrate(igrand,
          lower = 0, upper = 0.94574893543752026,
          y=y, mu=mu, phi=phi, p=p)$value


## PRE_Acceleration
area1 = integrate(igrand,
                  lower = 0.94574893543752026, upper = 9.8725536234497966     ,
                  y=y, mu=mu, phi=phi, p=p)$value

## TAIL
areaA = integrate(igrand,
                  lower = 9.8725536234497966     , upper = Inf,
                  y=y, mu=mu, phi=phi, p=p)$value

areaT = area0 + area1 + areaA 
cat("Area0 = ", area0, "\n")
cat("Area1 = ", area1, " (all done in one hit)\n")
cat("AreaA = ", areaA, " (all done in one hit)\n")
cat("=== TOTAL: ", areaT, " (all done in one hit)\n")



##################################

# library(tweedie); ptweedie.inversion(1, mu=1.4, phi=0.74, power=3, exact=TRUE)
source("./tweedieFunctions.R")

y=0.1 

## INITIAL
area0 = integrate(igrand,
                  lower = 0, upper = 19.604593198098087 ,
                  y=y, mu=mu, phi=phi, p=p)$value


## PRE_Acceleration
area1 = integrate(igrand,
                  lower = 19.604593198098087 , upper = 134.78996181244071      ,
                  y=y, mu=mu, phi=phi, p=p)$value

## TAIL
areaA = integrate(igrand,
                  lower = 134.78996181244071      , upper = Inf,
                  y=y, mu=mu, phi=phi, p=p)$value

areaT = area0 + area1 + areaA 
cat("Area0 = ", area0, "\n")
cat("Area1 = ", area1, " (all done in one hit)\n")
cat("AreaA = ", areaA, " (all done in one hit)\n")
cat("=== TOTAL: ", areaT, " (all done in one hit)\n")

areaT = area0 + area1 + areaA 
cat("Area0 = ", area0, "\n")
cat("Area1 = ", area1, " (all done in one hit)\n")
cat("AreaA = ", areaA, " (all done in one hit)\n")
cat("=== TOTAL: ", areaT, " (all done in one hit)\n")

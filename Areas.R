# library(tweedie); ptweedie.inversion(1, mu=1.4, phi=0.74, power=3, exact=TRUE)
source("./tweedieFunctions.R")

y=1 
mu=1.4
phi=0.74
p=3




mu <- 1
phi <- 4
p <- 1.5
# y <- 0.01
y <- 3




## INITIAL
area0 = integrate(igrand,
          lower = 0, upper = 1.1096230410971726,
          y=y, mu=mu, phi=phi, p=p)$value


## PRE_Acceleration
area1 = integrate(igrand,
                  lower = 1.1096230410971726, upper = 2.1314528304773006     ,
                  y=y, mu=mu, phi=phi, p=p)$value

## TAIL
areaA = integrate(igrand,
                  lower = 2.1314528304773006     , upper = 3.1672637245678148,
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

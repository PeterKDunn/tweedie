qtweedie <- function(p, xi = NULL, mu, phi, power = NULL){
  
  ### BEGIN preliminary work
  
  # SORT OUT THE NOTATION (i.e., xi VS power)
  if (verbose) cat("- Checking notation\n")
  out <- sort_Notation(xi = xi, power = power)
  xi <- out$xi
  power <- out$power
  xi.notation <- out$xi.notation
  index.par <- out$index.par
  index.par.long <- out$index.par.long ### MAY NOT BE NEEDED!!!

  
  # CHECK THE INPUTS ARE OK AND OF CORRECT LENGTHS
  if (verbose) cat("- Checking, resizing inputs\n")
  out <- check_Inputs(p, mu, phi, power)
  mu <- out$mu
  phi <- out$phi
  f <- array(0,
             dim = length(q) )
  if (details) regions <- array(0, dim = length(q))
  
  
  # IDENTIFY SPECIAL CASES
  if (verbose) cat("- Checking for special cases\n")
  out <- special_Cases(p, mu, phi, power)
  f <- out$f
  special_p_Cases <- out$special_p_Cases
  if (verbose & special_p_Cases) cat("  - Special case for p used\n")
  special_y_Cases <- out$special_y_Cases
  if (verbose & any(special_y_Cases)) cat("  - Special cases for first input found\n")
  
  ### END preliminary work
  
  
  
  
  
  len <- length(p) 
  
  # Some monkeying around to explicitly account for the cases p=1 and p=0
  ans <- ans2 <- rep( NA, length = len )
  if ( any(p == 1) ) ans2[p == 1] <- Inf
  if ( any(p == 0) ) ans2[p == 0] <- 0
  
  ans     <-  ans[ ( (p > 0) & (p < 1) ) ]
  mu.vec  <-  mu[ ( (p > 0) & (p < 1) ) ]
  phi.vec <-  phi[ ( (p > 0) & (p < 1) ) ]
  p.vec   <- p[ ( (p > 0) & (p < 1) ) ]
  
  for (i in (1 : length(ans)) ) {
    
    mu.1 <- mu.vec[i]
    phi.1 <- phi.vec[i]
    p.1 <- p.vec[i]  # This is the  qtweedie()  input p (a probability)
    pwr <- power[i]  # This is the Tweedie power, xi
    
    prob <- p.1 # Rename p to avoid confusion with  pwr: This is the  qtweedie()  input p (a probability)
    
    if ( pwr < 2 ) {
      qp <- qpois(prob, 
                  lambda = mu.1 / phi.1)
      if ( pwr == 1 ) ans[i] <- qp   
    }
    
    qg <- qgamma(prob,  
                 rate = 1 / (phi.1 * mu.1), 
                 shape = 1 / phi.1 )
    if ( pwr == 2 ) ans[i] <- qg
    
    # Starting values
    # - for 1<pwr<2, linearly interpolate between Poisson and gamma
    if ( (pwr > 1) & ( pwr < 2) ) {
      start <- (qg - qp) * pwr + (2 * qp - qg)
    }
    
    # - for pwr>2, start with gamma
    if ( pwr > 2 ) start <- qg
    
    # Solve!
    if ( ( pwr > 1) & (pwr < 2) ) { # This gives a *lower* bound on the value of the answer (if y>0)
      step <- dtweedie(y = 0, 
                       mu = mu.1, 
                       phi = phi.1, 
                       power = pwr)
      # This is P(Y = 0), the discrete "step"
      
      if ( prob <= step ) {
        ans[i] <- 0
      }
    }
    
    if ( is.na(ans[i]) ) { # All cases except Y=0 when 1 < pwr < 2
      
      
      pt2 <- function( q, 
                       mu, 
                       phi, 
                       pwr, 
                       p.given = prob ){ 
        
        ptweedie(q = q, 
                 mu = mu, 
                 phi = phi, 
                 power = pwr ) - p.given
      }
      
      pt <- pt2( q = start, 
                 mu = mu.1, 
                 phi = phi.1, 
                 pwr = pwr,
                 p.given = prob)
      
      #      cat("*** Start   =",start,"; pt =",pt,"\n")
      
      if ( pt == 0 ) ans2[i] <- start
      
      if ( pt > 0 ) { 
        loop <- TRUE
        start.2 <- start
        #start.2 <- 0   # Perhaps set this if too many attempts otherwise
        while ( loop ) {
          # Try harder
          start.2 <- 0.5 * start.2
          if (pt2( q = start.2, 
                   mu.1, 
                   phi.1, 
                   pwr, 
                   p.given = prob )<0 ) loop = FALSE
          #      cat(">>> Start.2 =",start.2,"; pt = ",pt2( q=start.2, mu.1, phi.1, pwr, p.given=prob ),"\n")
          # RECALL:  We are only is this part of the loop if  pt>0
        }
      }
      
      #      cat("*** Start.2 =",start.2,"\n")
      if ( pt < 0) {
        loop <- TRUE
        start.2 <- start
        
        while ( loop ) {
          # Try harder
          start.2 <- 1.5 * (start.2 + 2)
          if (pt2( q = start.2, 
                   mu.1, 
                   phi.1, 
                   pwr, 
                   p.given = prob ) > 0 ) loop = FALSE
          # RECALL:  We are only is this part of the loop if  pt<0
        }
      }
      
      #cat("start, start.2 =",start, start.2,"\n")
      #cat("pt2(start, start.2) =",
      #   pt2(start, mu=mu.1, phi=phi.1, pwr=pwr, p.given=prob), 
      #   pt2(start.2, mu=mu.1, phi=phi.1, pwr=pwr, p.given=prob),"\n")
      
      out <- uniroot(pt2, 
                     c(start, start.2), 
                     mu = mu.1, 
                     phi = phi.1, 
                     p = pwr, 
                     p.given = prob )
      #print(out)
      
      ans[i] <- uniroot(pt2, 
                        c(start, start.2), 
                        mu = mu.1, 
                        phi = phi.1, 
                        p = pwr, 
                        p.given = prob, 
                        tol = 0.000000000001 )$root
    }
    
  }
  
  ans2[ is.na(ans2) ] <-  ans
  ans2
}



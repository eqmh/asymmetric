#' covstop: a function to implement coverage-based stopping
#' 
#' @param comm a site (rows)-by-species (columns) data.frame
#' @param type whether the output should be in number of individuals ("individuals") or in
#' number of sampling units ("samples")
#' @param Cn desired level of coverage
#' @param se whether bootstrapped standard errors should be returned
#' @param knots a number of knots of computation, default is 40
#' 
#' @return the number of individuals or samples required to meet the desired level of coverage
#' 
#' @examples
#' comm <- sample(100)
#' 
#' comm[sample(1:100, 50, replace = F)] <- 0
#' 
#' comm <- matrix(comm, nrow = 10)
#' 
#' covstop(comm)
#' 
covstop <- function(comm, type = "samples", Cn = 0.90, se = FALSE, knots = 40) {
  
  if(Cn > 1) stop("Coverage level must be between (0, 1)")
  
  if(type == "samples") comm <- data.frame(N = nrow(comm), t(colSums(comm > 0)))
  
  comm <- split(comm, 1:nrow(comm))
  
  sc <- sapply(comm, function(Spec) {
    
    Spec <- as.numeric(Spec)
    
    endpoint <- 2 * sum(Spec)
  
    if(type == "individuals") {
      
      n <- sum(Spec)
      
      if (endpoint <= n) {
        
          m <- floor(seq(1, endpoint, length = floor(knots)))
          
        } else {
          
          m <- c(floor(seq(1, sum(Spec) - 1, length.out = floor(knots/2) - 1)), sum(Spec), floor(seq(sum(Spec) + 1, to = endpoint, length.out = floor(knots/2))))
          
          }
      
      m <- c(1, m[-1])
      
      C <- Chat.Ind(Spec, m)
      
      # if(se == TRUE) {
      #   
      #   Prob.hat <- EstiBootComm.Ind(Spec)
      #   
      #   Abun.Mat <- rmultinom(nboot, n, Prob.hat)
      #   
      #   error.C <- qnorm(1 - (1 - 0.95) / 2) * apply(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)), 1, sd, na.rm = TRUE)
      #   
      # }
        
    } else if(type == "samples") {
      
      nT <- Spec[1]
      
      if(endpoint <= nT) {
        
          m <- floor(seq(1, endpoint, length.out = floor(knots)))
          
          } else {
            
            m <- c(floor(seq(1, nT - 1, length.out = floor(knots/2) - 1)), nT, floor(seq(nT + 1, to = endpoint, length.out = floor(knots/2))))
            
            }
      
        m <- c(1, m[-1])
 
      C <- Chat.Sam(Spec, m)
      
    }
    
    
    m[which.min(C <= Cn)]
    
  } )
  
  # if(se == TRUE) 
  #   
  #   paste0("Number of ", type, " needed to reach ", Cn*100, "% coverage: ", as.numeric(sc[which.max(sc)]), "+/-", as.numeric(error.C))
  
  paste0("Number of ", type, " needed to reach ", Cn*100, "% coverage: ", as.numeric(sc[which.max(sc)]))
  
}

### helper functions from iNEXT 
### credit: T.C. Hsieh
Chat.Ind <- function (x, m) {
  x <- x[x > 0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 
                                                            1)/n * f1^2/2/f2)
  A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
  Sub <- function(m) {
    if (m < n) {
      xx <- x[(n - x) >= m]
      out <- 1 - sum(xx/n * exp(lgamma(n - xx + 1) - lgamma(n - 
                                                              xx - m + 1) - lgamma(n) + lgamma(n - m)))
    }
    if (m == n) 
      out <- 1 - f1/n * A
    if (m > n) 
      out <- 1 - f1/n * A^(m - n + 1)
    out
  }
  sapply(m, Sub)
}

Chat.Sam <- function(x, t){
  nT <- x[1]
  y <- x[-1]
  y <- y[y>0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)  #estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  Sub <- function(t){
    if(t < nT) {
      yy <- y[(nT-y)>=t]
      out <- 1 - sum(yy / U * exp(lgamma(nT-yy+1)-lgamma(nT-yy-t+1)-lgamma(nT)+lgamma(nT-t)))     
    }
    #if(t < nT) out <- 1 - sum(y / U * exp(lchoose(nT - y, t) - lchoose(nT - 1, t)))
    if(t == nT) out <- 1 - Q1 / U * A
    if(t > nT) out <- 1 - Q1 / U * A^(t - nT + 1)
    out
  }
  sapply(t, Sub)		
}

EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0)   #observed species
  n <- sum(Spec)        #sample size
  f1 <- sum(Spec == 1)   #singleton 
  f2 <- sum(Spec == 2)   #doubleton
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  a <- f1/n*A
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  if(f0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(f0.hat), ceiling(f0.hat))  	#estimation of relative abundance of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))		  							#Output: a vector of estimated relative abundance
}

EstiBootComm.Sam <- function(Spec)
{
  nT <- Spec[1]
  Spec <- Spec[-1]
  Sobs <- sum(Spec > 0)   #observed species
  Q1 <- sum(Spec == 1) 	#singleton 
  Q2 <- sum(Spec == 2) 	#doubleton
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)	#estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  a <- Q1/nT*A
  b <- sum(Spec / nT * (1 - Spec / nT) ^ nT)
  
  if(Q0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  
  Prob.hat <- Spec / nT * (1 - w * (1 - Spec / nT) ^ nT)					#estimation of detection probability of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(Q0.hat), ceiling(Q0.hat))  	#estimation of detection probability of unseen species in the sample
  return(sort(c(Prob.hat, Prob.hat.Unse), decreasing=TRUE))									#Output: a vector of estimated detection probability
}

#' estimateD.new = faster coverage-based rarefaction 
#' @param x a community (rows) by species (cols) matrix
#' @param q the level of q to be used for calculating diversity (default is richness)
#' @param level a user-defined level of coverage
#' @return a data.frame with the sample size, method of estimation, sample coverage, and estimated diversity
#' 
estimateD.new <- function(x, q = 0, level=NULL){
  
  C <- level
  
  if(is.null(C))
    
    C <- min(unlist(apply(x, 1, function(x) Chat.Ind(x, sum(x)))))
  
  do.call(rbind, apply(x, 1, function(x) invChat.Ind(x, q, C)))
  
}

#' helper functions
Chat.Ind <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)		
}

invChat.Ind <- function(x, q, C){
  
  m <- NULL # no visible binding for global variable 'm'
  
  n <- sum(x)
  
  refC <- Chat.Ind(x,n)
  
  f <- function(m, C) abs(Chat.Ind(x,m)-C)
  
  if(refC > C) {
    opt <- optimize(f, C=C, lower=0, upper=sum(x))
    mm <- opt$minimum
    mm <- round(mm)
  }
  
  if(refC <= C) {
    f1 <- sum(x==1)
    f2 <- sum(x==2)
    if(f1>0 & f2>0){A <- (n-1)*f1/((n-1)*f1+2*f2)}
    if(f1>1 & f2==0){A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
    if(f1==1 & f2==0){A <- 1}
    if(f1==0 & f2==0){A <- 1}
    mm <- (log(n/f1)+log(1-C))/log(A)-1
    mm <- n+mm
    mm <- round(mm)
  }
  
  if(mm > 2*n) warning("The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias.")
  
  method <- ifelse(mm<n, "interpolated", ifelse(mm==n, "observed", "extrapolated"))
  
  out <- data.frame(m=mm,
                    method=method, 
                    coverage=round(Chat.Ind(x,mm),3),
                    diversity=round(Dqhat.Ind(x,q,mm),3)
  )
  
  out
  
}

Dqhat.Ind <- function(x, q, m){
  x <- x[x > 0]
  n <- sum(x)
  
  fk.hat <- function(x, m){
    x <- x[x > 0]
    n <- sum(x)
    if(m <= n){
      Sub <- function(k)	sum(exp(lchoose(x, k) + lchoose(n - x, m -k) - lchoose(n, m)))
      sapply(1:m, Sub)
    }
    
    else {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      A <- ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
      C.hat <- 1 - f1 / n * A
      p.hat <- x / n * C.hat			
      Sub <- function(k)	sum((choose(m, k) * p.hat^k * (1 - p.hat)^(m - k)) / (1 - (1 - p.hat)^n))
      sapply(1:m, Sub)
    }
  }
  
  D0.hat <- function(x, m){
    x <- x[x > 0]
    n <- sum(x)
    Sub <- function(m){
      if(m <= n){
        Fun <- function(x){
          if(x <= (n - m)) exp(lgamma(n - x + 1) + lgamma(n - m + 1) - lgamma(n - x - m + 1) - lgamma(n + 1))
          else 0
        }
        sum(1 - sapply(x, Fun))
      }
      else {
        Sobs <- sum(x > 0)
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)	#estimation of unseen species via Chao1
        A <- n*f0.hat/(n*f0.hat+f1)
        ifelse(f1 ==0, Sobs ,Sobs + f0.hat * (1 - A ^ (m - n)))	
      }
    }
    sapply(m, Sub)
  }
  
  D1.hat <- function(x, m){
    x <- x[x > 0]
    n <- sum(x)
    Sub <- function(m){
      if(m < n){
        k <- 1:m
        exp(-sum(k / m * log(k / m) * fk.hat(x, m)))
      }
      else{
        #UE=sum(sapply(1:(n-1),function(k){sum(1/k*x/n*exp(lchoose(n-x,k)-lchoose(n-1,k)))}))
        UE <- sum(x/n*(digamma(n)-digamma(x)))
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        A <- 1 - ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
        #A=2*sum(x==2)/((n-1)*sum(x==1)+2*sum(x==2))
        B <- ifelse(A<1, sum(x==1)/n*(1-A)^(-n+1)*(-log(A)-sum(sapply(1:(n-1),function(k){1/k*(1-A)^k}))), 0)
        D.hat <- exp(UE+B)
        Dn <- exp(-sum(x / n * log(x / n)))
        
        a <- 1:(n-1)
        b <- 1:(n-2)
        Da <-  exp(-sum(a / (n-1) * log(a / (n-1)) * fk.hat(x, (n-1))))
        Db <-  exp(-sum(b / (n-2) * log(b / (n-2)) * fk.hat(x, (n-2))))
        Dn1 <- ifelse(Da!=Db, Dn + (Dn-Da)^2/(Da-Db), Dn)
        b <- ifelse(D.hat>Dn, (Dn1-Dn)/(D.hat-Dn), 0)
        # b <- A
        ifelse(b!=0, Dn + (D.hat-Dn)*(1-(1-b)^(m-n)), Dn)
      }
    }
    sapply(m, Sub)
  }
  
  D2.hat <- function(x, m){
    Sub <- function(m) 1 / (1 / m + (1 - 1 / m) * sum(x * (x - 1) / n / (n - 1)))
    sapply(m, Sub)
  }
  
  Dq.hat <- function(x, m){
    Sub <- function(m){
      k <- 1:m
      sum( (k / m)^q * fk.hat(x, m))^(1 / (1 - q))
    }
    sapply(m, Sub)
  }
  if(q == 0) D0.hat(x, m)
  else if(q == 1) D1.hat(x, m)
  else if(q == 2) D2.hat(x, m)
  else Dq.hat(x, m)
}
#' Covstop: a function to implement coverage-based stopping
#' 
#' @param comm a species-by-site data.frame or list
#' @param Cn desired level of coverage
#' 
#' @return the sample size (number of individuals) required to meet the desired
#' level of coverage
#' 
covstop <- function(comm, Cn = 0.90) {
  
  if(Cn > 1) stop("Coverage level must be between (0, 1)")
  
  if(!is.list(comm)) comm <- split(comm, 1:ncol(comm))
  
  sapply(comm, function(Spec) {
    
    n <- sum(Spec)
    
    endpoint <- knots <- 2*sum(Spec)
    
    # if (is.null(m)) {
      if (endpoint <= n) {
        m <- floor(seq(1, endpoint, length = floor(knots)))
      } else {
        m <- c(floor(seq(1, sum(Spec) - 1, length.out = floor(knots/2) - 
                           1)), sum(Spec), floor(seq(sum(Spec) + 1, to = endpoint, 
                                                     length.out = floor(knots/2))))
      }
      m <- c(1, m[-1])
    # }
    
    Chat.Ind <- function (x, m) 
    {
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
    
    cc <- cbind.data.frame(
      n = m,
      C = Chat.Ind(Spec, m)
    )
    
    cc[which.min(cc$C <= Cn), "n"]
    
  } )
  
}

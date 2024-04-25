require(condMVNorm)
get_fraction_m_alt <- function(meanA=0.5, m=2, n=3, ...){
  #n x n matrix
  A   <- diag(n) + meanA
  A   <- set_diagonal(A, 1)
  Am <- A %*% diag(c(rep(1,m),rep(0,n-m)))
  #Am2<- diag(c(rep(1,m),rep(0,n-m))) %*% Am
  diag(Am) <- c(rep(1,m), rep(-1,n-m)) 
  #diag(Am2)<- 1
  B        <- diag(c(rep(1, n))) #constraints
  #feas1 <- 2^(n)*calculate_omega_constraint(A=Am, B=B)^(n)
  feas1 <- 2^(n)*calculate_omega(Am)^(n)
  print(feas1)
  Ainv <- solve(A)
  # Calculate and partition mu and sigma
  mu <- c(Ainv%*%rep(0, n))
  sigma <- Ainv%*%diag(rep(1^2, n))%*%t(Ainv)
  mu1 <- mu[1:m]; mu2 = mu[(m+1):n]
  sigma11 <- sigma[c(1:m), c(1:m)]
  sigma12 <- sigma[c(1:m), c((m+1):n)]
  sigma21 <- sigma[c((m+1):n), c(1:m)]
  sigma22 <- sigma[c((m+1):n), c((m+1):n)]
  muBar <- mu1 + sigma12 %*% solve(sigma22) %*% (0-mu2)
  sigmaBar <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  feas2 <- 2^n*pmvnorm(lower=rep(0,m), upper=Inf, mean=as.vector(muBar), 
                   sigma=sigmaBar)
  #feas2 <- (2^n)*pmvnorm(lower=rep(0,n), upper=Inf, mean=as.vector(mu), 
  #                 sigma=sigma)
  print(feas2)
}






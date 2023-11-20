library(deSolve)
library(mvtnorm)
library(tidyverse)
library(igraph)
library(pracma)
library(gridExtra)
library(feasoverlap)
library(reticulate) #to run python code from R
use_python("/usr/local/bin/python3.10") #load your preferred python version
source_python("MVN.py")
extinctionThreshold <- 1e-3

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#

# Functions ------------------------------------------------------------
## To make all sorts of matrices -------------------------------------
#Make block diagonal matrix by stacking p times A
# on the diagonal. The matrix A can be simply perfectly copied (vary=0),
# or varied to some extend among copies, as quantified by  
# "vary" = the width of the uniform (max-min) around the elements of A
make_block_diagonal <- function(A, p, vary = 0, ...) {
  A_temp <- array(0, dim=dim(A)*p) # all zeros
  for (i in 1:p) {
    A_vary <- array(data = runif(prod(dim(A)), 1 - vary, 1 - vary), dim = dim(A)) #local deviation of A in block i
    A_temp[(ncol(A)*i-ncol(A)+1):(ncol(A)*i), (ncol(A)*i-ncol(A)+1):(ncol(A)*i)] <- A * A_vary
    }
  return(A_temp)
}

set_diagonal <- function(A, d){
  diag(A) <- d
  A
}

#Make a dispersal (D) matrix
# all diagonal blocks are zero
# off diagonals blocks are diagonal matrices containing the sp-specific dispersal rates
# if two sites are connected, they exchange individuals in both directions (but allowing for different rates)
# d = dispersal rate 
# p = nr of locations
# n = nr of sp
# connectivity = *fraction* of all possible combos that are connected
make_D <- function(p=3, n=2, connectivity=1, d=1e-3, ...){
  D    <- array(0, dim=c(n*p, n*p))#matrix w correct dimensions and all zeros
  pos_links<- combn(c(1:p), 2) # all possible combos
  nr_links <- round(connectivity * ncol(pos_links), 0) #nr of links between locations
  links    <- pos_links[,sample(c(1:ncol(pos_links)), nr_links)]#sample that nr of links
  for (i in 1:ncol(links)) {
    D_12 <- diag(rep(d, n))#make random dispersal matrix from 1 to 2
    D_21 <- diag(rep(d, n))#make random dispersal matrix from 2 to 1
    pos_1 <- (n*links[1,i]-n+1):(n*links[1,i]) #position of first location
    pos_2 <- (n*links[2,i]-n+1):(n*links[2,i]) #position of second location
    D[pos_1, pos_2]<-D_12
    D[pos_2, pos_1]<-D_21
  }
  return(D)
}

#make symmetric and ditch diagonal
make_symmetric <- function(A){
  A*upper.tri(A) + t(A)*lower.tri(A)
}

#generates p intrinsic growth rate vectors, stacked underneath each other
#for use in spatial LV simulations
#negative_sign is a vector telling which species needs to have a negative intrinsic growth rate
#example: c(2,3) is the case where the intrinsic growth rate of sp 2 and 3 carry a negative sign
#If you want none to have a negative r, just go negative_sign= some i>n
#n is nr of species, p is nr of patches
make_R_spatial <- function(n, p, negative_sign=10000, ...){
  #sample from a sphere to get the directions right
  rawRs <- t(sample_sphere_surface(dim=n, n = p, radius = 1, positive = TRUE))
  #now make sure the mean r across species is 1 at all patches
  rawRs <- rawRs %*% diag(1/colMeans(rawRs))
  rawRs <- diag(1/rowMeans(rawRs)) %*% rawRs  
  #and put everything into a nice format
  Rs     <- as_tibble(rawRs)%>%
    rename_with(~gsub("V","", .x, fixed = TRUE)) %>% #give correct sign to intrinsic growth rate of consumer
    mutate(location=c(1:p)) %>%
    pivot_longer(!location) %>%
    mutate(value=value*(-1)*name%in%negative_sign + value*(1-name%in%negative_sign))
  return(Rs$value)
}

## LV -------------------------------------

#Cut-off function to limit the values of the state to enhance 
#numerical stability when solving the odes
cutoff <- function(n, a) {
  return(1*(n<a)*(3*(n/a)^2-2*(n/a)^3)*ifelse(n<0, 0, 1)+ifelse(n<a, 0, 1))
}

#the LV model; R= a vector of max. intrinsic growth rates; 
#A= matrix with linear species interactions, global var 
run_LV <- function(time, state, pars) {
  return(list(as.numeric(
    ((state*(pars$R+pars$A%*%state)))*cutoff(state, 1e-10)
  )))
}

#Spatial LV model; R= a vector of max. intrinsic growth rates; 
#A= matrix with linear species interactions, D=dispersal matrix
run_LV_spatial <- function(time, state, pars) {
  return(list(as.numeric(
    ((state*((pars$R-colSums(pars$D))-pars$A%*%state))+pars$D%*%state)*
      cutoff(state, 1e-10)
  )))
}
#compute NHat from spatial LV simulation
#n = nr of sp
#R, A, D: as in run_LV_spatial
#max_time = max simulation time

get_NHat <- function(n, R, A, D, max_time=100, N0, ...){
  p         <- length(R)/n #infer nr of patches
  densities <- ode(func=run_LV_spatial, y=N0, 
                   parms=list(R=R, A=A, D=D), 
                   times=seq(1,max_time,0.1))
  locations <- NULL; for (i in 1:p) {locations <- c(locations, rep(i, n))}
  return(as_tibble(cbind(density=tail(densities, 1)[-1], 
                         location=locations,
                         sp=rep(c(1:n), p)))) 
}

## Others ------
get_N_total <- function(meanA=0.2, d=1, n=10, r=1){ #a=comp.strength; d=self limit.
  n*r/(d+meanA*(n-1))
}

#compute fraction of patches with m species in a metacommunity of n species
#without dispersal. For m>1
get_fraction_m <- function(meanA=0.5, m=2, n=3, ...){
  #n x n matrix
  A   <- diag(n) + meanA
  A   <- set_diagonal(A, 1)
  Am <- A %*% diag(c(rep(1,m),rep(0,n-m)))
  #Am2<- diag(c(rep(1,m),rep(0,n-m))) %*% Am
  diag(Am) <- c(rep(1,m), rep(-1,n-m)) 
  #diag(Am2)<- 1
  B        <- diag(c(rep(1, n))) #constraints
  feas1combo <- 2^(n)*calculate_omega_constraint(A=Am, B=B)^(n)
  choose(n, m) * feas1combo
  }

#get Ni when N0i is larger than zero for all i
get_Ni_N0iLargerThen0 <- function(a=0.5, n=10, d=1e-4, p=100, N0Inv=1/1, ri=1, NTotalK=10) {
  (1/(1 + a*(-1 + n)))*((a + a*n*(-1 + ri) + ri - 2*a*ri)/(1 - a) + 
   d*(1 - p + ((-1 + a*(4 - 2*n + a*(-5 + a*(-2 + n)*(-1 + n) - N0Inv - n*(-5 + n + (-2 + n)*N0Inv)) + (1 + a*(-2 + n))*(-1 + n)*N0Inv*ri))*NTotalK)/
        ((-1 + a)*(a - a*n + ri + a*(-2 + n)*ri))))
}

#Simplified version, assuming that N0i times the mean of 1/N0i across species is roughly equal to 1
#N0i has no default because needs to be computed by get_N0i
get_Ni_N0iLargerThen0 <- function(a=0.5, n=10, d=1e-4, p=100, 
                                  N0i, NTotalK=10) {
  ((1 + a*(-1 + n))*N0i^2 + d*(N0i - N0i*p + NTotalK))/(N0i + a*(-1 + n)*N0i)
}

#get the density of a species i when there is no dispersal
get_N0i <- function(a=0.5, n=10, r=1, ri=0.1){
  (a*(n - 1)*r - ri*(a*(n - 2) + 1))/((a - 1)*(a*(n - 1) + 1))
}

#get Ni when N0i is equal to zero (for at least the focal sp i)
#SumN0j is a RV: it is the sum of biomass across species in a patch, 
#where m species (m<n) coexist when there is no dispersal.
get_Ni_N0iEqualTo0 <- function(d=1e-4, NTotalK=10, ri=1, meanA=0.5, SumN0j){
  (d*NTotalK)/(ri - meanA*SumN0j)
}

#make a distribution of nr of species in a patch in
#case there is no dispersal
make_distribution <- function(n, meanA=0.5){
  probs <- NULL
  for (m in c(1:n))
  {
    probs <- c(probs, get_fraction_m(meanA=meanA, m=m, n=n))
  }
  probs
}
#truncate a distribution (i.e. cut off fraction below (if ditch="down") 
#or above (ditch="up") a certain quantile q)
#pdfFitted = output of density()
#q = quantile, where 1 = 100th quantile, 0.5 = 50th etc..
trunc_dist <- function(pdfFitted, q=0.5, ditch="up"){
  cumprob <- cumsum(pdfFitted$y)/sum(pdfFitted$y) #discretize to cum. proba.
  #identify x to be ditched
  if(ditch=="down") {Qs<-which(cumprob>=q)} else {Qs<-which(cumprob<=q)}
  list(x=pdfFitted$x[Qs], y=pdfFitted$y[Qs])
}

#get mean of a truncated pdf
#pdfFitted = output of trunc_dist()
#q = quantile, where 1 = 100th quantile, 0.5 = 50th etc..
get_mean_trunc <- function(pdfFitted, q=0.5, ditch="up"){
  pdfTrunc <- trunc_dist(pdfFitted, q=q, ditch=ditch)
  sum(pdfTrunc$x*pdfTrunc$y)/sum(pdfTrunc$y)
}

#mean R of the persisting sp
#x = factor to divide a*Nt by in order to get the mean
get_RMeanM <- function(a=0.8, m=1, n=4, r=1.01, x=2.2){
  ((1 + a*(-1 + m))*n*r*x)/(m*(a*(n + m*(-1 + x) - x) + x)) 
} 
#mean R of the persisting sp, complicated version
#RMeanM <- function(a=0.8, m=1, n=4, r=1.01){
#  ((sqrt(3)*sqrt((m + a*(-1 + m)*m)^2*(3*(2 + a*(-2 + m + n))^2 + 2*a*n*(6 + a*(-6 + m + 5*n))*r + 3*a^2*n^2*r^2)) + 
#       3*(1 + a*(-1 + m))*m*(2 + a*(-2 + m + n - n*r)))/(2*a*m^2*(3 + a*(-3 + m + 2*n))))
#} 
#probability of extinction w/o dispersal (for visualisation purposes)
#feas = vector w feasibilities of all m for which 1<=m<=n
get_prob_ext_wo <- function(n, feas){
  n <- length(feas)
  probs <- tibble(n=n, m=c(1:n), feas=feas) %>%
    mutate(prob = feas*(1-m/n))
  sum(probs$prob)
}

#value of N1i when persistence w/o dispersal
# N0i = density of i w/o dispersal
# m = nr of species coexisting w/o dispersal
# n = nr of species regionally available
# N1 = mean of N1j across j that don't persist
# p = nr of patches
# SumN0k = total density across whole network for a species (same for all sp)
# N0Inv = mean of the inverse of N0i
get_N1i_no_ext_wo <- function(meanA, N0i, m, n, N1, p, SumN0k, 
                              N0Inv){
  (-((-1 + meanA)*N0i*(1 + meanA*(m - n)*N1 - p)) + SumN0k + 
     meanA*(-2 + m + N0i*N0Inv - m*N0i*N0Inv)*SumN0k)/((-1 + meanA)*(1 + meanA*(-1 + m))*N0i)
}




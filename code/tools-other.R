library(mvtnorm)
library(pracma)
library(gridExtra)
library(feasoverlap)
library(reticulate) #to run python code from R
use_python("/usr/local/bin/python3.10") #load your preferred python version
source_python("MVN.py")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#

get_N_total <- function(meanA=0.2, d=1, n=10, r=1){ #a=comp.strength; d=self limit.
  n*r/(d+meanA*(n-1))
}

#get the density of a species i when there is no dispersal
get_N0i <- function(a=0.5, n=10, r=1, ri=0.1){
  (a*(n - 1)*r - ri*(a*(n - 2) + 1))/((a - 1)*(a*(n - 1) + 1))
}

#get N1i when i is excluded w/o dispersal
#NTotalK = total abundance of i across all patches without the focal patch
#ri = r of i, given that it gets excluded w/o dispersal
#a = interaction strength
#SumN0j = abundance of all species 
get_N1iExc <- function(NTotalK=10, ri=1, a=0.5, SumN0j){
  NTotalK/(a*SumN0j - ri)
}

#get N1i when i persists w/o dispersal
#a = interaction strength
#n and m: nr of species in the regional pool and in the focal patch
#meanN1iExc = mean of N1iExc across species
#rho = mean of rhoi across species
#rhoi = rho for species i, given that it persists w/o dispersal
get_N1iPer <- function(a=0.5, n=4, m=1, rho, rhoi, p, meanN1iExc){
  (a*((1-a)*meanN1iExc*(n-m)+(m-1)*rho+(2-m)*rhoi) + (p-1)*(1-a) - rhoi)/
    ((1-a)*(a*(m-1)+1))
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

#sample ri's from a truncated dist
# PDF = a pdf, with $x the values and $y the density
# cutoff = threshold for truncation
# ditch = whether to ditch everything below or above
sample_ri <- function(samplesize=1, PDF, cutoff, ditch="below", ...){
  ditchBelow <- sample(x = PDF$x[which(PDF$x>cutoff)], size=1, 
                       prob = PDF$y[which(PDF$x>cutoff)], replace = T)
  ditchAbove <- sample(x = PDF$x[which(PDF$x<cutoff)], size=1, 
                       prob = PDF$y[which(PDF$x<cutoff)], replace = T)
  (ditch=="below")*ditchBelow + (ditch=="above")*ditchAbove
}

# Get patch occupancy
# Input:
# -fm: distribution of number of coexisting species across network
# -pPersist: probability for a species to persist if it does not persist without dispersal
# -n: total number (i.e. number of regionally available) species
# -Xi: feasibility for n species
get_patch_occupancy <- function(fm, pPersist, n, Xi, ...){
  sum(fm*pPersist^(n-c(1:n)))*(1-Xi) + Xi
}
 
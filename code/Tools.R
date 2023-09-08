library(deSolve)
library(mvtnorm)
library(tidyverse)
library(igraph)
library(pracma)
library(gridExtra)
library(reticulate) #to run python code from R
use_python("/usr/local/bin/python3.10") #load your preferred python version
source_python("MVN.py")

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
  rawRs <- diag(1/rowMeans(rawRs)) %*% rawRs 
  #abd put everything into a nice format
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
get_N_total <- function(mean_a=0.2, d=1, l=10, r=1){ #a=comp.strength; d=self limit.
  l*r/(d+mean_a*(l-1))
}






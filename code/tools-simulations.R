library(deSolve)
library(tidyverse)
library(igraph)
extinctionThreshold <- 1e-3

# Functions ------------------------------------------------------------
## To make all sorts of matrices -------------------------------------
#Make block diagonal matrix by stacking p times A
# on the diagonal. The matrix A can be simply perfectly copied (vary=0),
# or varied to some extend among copies, as quantified by  
# "vary" = the width of the uniform (max-min) around the elements of A
make_block_diagonal <- function(A, p, vary = 0, ...) {
  A_temp <- array(0, dim=dim(A)*p) # all zeros
  for (i in 1:p) {
    A_vary <- array(data = runif(prod(dim(A)), 1 - vary, 1 + vary), dim = dim(A)) #local deviation of A in block i
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
#n is nr of species, p is nr of patches, k is fold increase of the strongest species' growth rate
make_R_spatial <- function(n, p, k=1, negative_sign=10000, ...){
  #sample from a sphere to get the directions right
  rawRs <- t(sample_sphere_surface(dim=n, n = p, radius = 1, positive = TRUE))
  #now make sure the mean r across species is 1 at all patches
  rawRs <- rawRs %*% diag(1/colMeans(rawRs))
  rawRs <- diag(1/rowMeans(rawRs)) %*% rawRs  
  rawRs <- rawRs %*% diag(c(k, rep(1/k, n-1))) # relax regional equivalence: make one species more competitive than expected across all patches
  colnames(rawRs) <- c(1:n)
  #and put everything into a nice format
  Rs     <- as_tibble(rawRs)%>%
    mutate(location=c(1:p)) %>%
    pivot_longer(!location) %>%
    mutate(value=value*(-1)*name%in%negative_sign + value*(1-name%in%negative_sign))
  return(Rs$value)
}

# Generate random landscape with randomly drawn patch coordinates, each between 0 and 1.
# Input:
# - nPatch: Number of patches in the landscape.
# - dim: Landscape dimensionality (1 for a lakeshore; 2 for an area, etc.).
# Output:
# - An nPatch x dim data frame whose (i,j)th entry is the jth coordinate of patch i.
#   The rows are unnamed; the columns are named "V1", "V2", ..., "V<dim>".
make_randomCoords <- function(nPatch, dim = 2) {
  matrix(runif(nPatch * dim, 0, 1), nrow = nPatch, ncol = dim) |>
    as.data.frame(row.names = NULL)
}

# Create dispersal matrix.
# Input:
# - coords: Patch coordinates (either as a matrix, data frame, or tibble),
#           with patches in the rows and coordinate values in the columns.
# - kernel: A function specifying the functional form of the dispersal kernel. The
#           function should have a single argument (distance between patches).
# Output:
# - A matrix, with diagonal entries zeroed out, where the (i,j)th entry is the
#   dispersal rate from patch j to patch i.
dispMatrix <- function(coords, kernel) {
  coords |> # Start with table of patch coordinates
    dist() |> # Compute pairwise distances between all possible combinations
    as.matrix() |> # Convert result from a "dist" object to a regular matrix
    kernel() |> # Apply the dispersal kernel function to each entry of the matrix
    unname() |> # Remove unnecessary row and column names from the matrix
    (`diag<-`)(0) # Set all diagonal entries to 0
}


# Create a list of dispersal matrices over the same landscape. Each matrix is generated
# with a specified dispersal kernel function.
# Input:
# - coords: Patch coordinates (either as a matrix, data frame, or tibble),
#           with patches in the rows and coordinate values in the columns.
# - kernelList: A list of functions specifying the functional forms of the dispersal
#               kernel, for each species. The functions should all have a single
#               argument (distance between patches).
# Output:
# - A list of dispersal matrices, each nrow(coords) x nrow(coords) in dimension. The
#   first list entry is the matrix for species 1, the second for species 2, and so on.
dispMatrixList <- function(coords, kernelList) {
  # For each function in kernelList, create a matrix; store them in a list:
  Map(\(f) dispMatrix(coords, f), kernelList)
}


# Create a matrix whose every entry is zero, except one single entry along its
# diagonal which is 1. This will be useful for creating a community-wide dispersal
# matrix using sums of Kronecker products.
# Input:
# - dim: The dimension of the (square) matrix.
# - n: The index of the entry along the diagonal that should be 1 instead of 0.
# Output:
# - A dim x dim matrix which is 1 at its (n,n)th entry and 0 otherwise.
oneHotMatrix <- function(dim, n) {
  rep(0, dim) |> # Vector of `dim` number of 0s
    replace(n, 1) |> # Set the nth entry to 1
    diag() # Create matrix with the above vector in the diagonal and zeros elsewhere
}


# Create S*L x S*L community-wide dispersal matrix, where S is the number of species
# and L is the number of patches.
# Input:
# - coords: Patch coordinates (either as a matrix, data frame, or tibble),
#           with patches in the rows and coordinate values in the columns.
# - kernelList: A list of functions specifying the functional forms of the dispersal
#               kernel, for each species. The functions should all have a single
#               argument (distance between patches).
# Output:
# - An S*L x S*L matrix, where S is the number of species and L is the number of patches.
dispMatrixCommunity <- function(coords, kernelList) {
  S <- length(kernelList) # Number of species
  L <- nrow(coords) # Number of patches
  dispMatrices <- coords |> dispMatrixList(kernelList) # S matrices in a list, each LxL
  fullDisp <- matrix(0, S*L, S*L) # Initialize full dispersal matrix
  # Add each species' contribution to the expanded dispersal matrix:
  for (i in 1:S) fullDisp <- fullDisp + dispMatrices[[i]] %x% oneHotMatrix(S, i)
  fullDisp # Return the assembled result
}

# Convenience function that works like `dispMatrixCommunity`, except it assumes that
# The dispersal kernels are all exponential. So instead of a list of dispersal functions,
# it receives a list of numbers - the characteristic dispersal distance per unit time
# for each species.
# Input:
# - coords: Patch coordinates (either as a matrix, data frame, or tibble),
#           with patches in the rows and coordinate values in the columns.
# - dispDistanceVector: Vector of dispersal distances. Its ith entry is for species i.
#                       Using these values, the kernels will be \(x) exp(x / d), where
#                       d is the appropriate entry of dispDistanceVector.
dispMatrixCommunityExp <- function(coords, dispDistanceVector) {
  # Create list of functions; e.g., list(\(x) exp(-x/1), \(x) exp(-x/2), \(x) exp(-x/3))
  # in case dispDistanceVector is given as c(1, 2, 3):
  Map(\(d) \(x) exp(-x / d), dispDistanceVector) |>
    # Now apply `dispMatrixCommunity` to these functions and the supplied coordinates:
    dispMatrixCommunity(coords = coords)
}

# Function to rescale the created D matrix to fix the mean 
# dispersal rate. Only makes sense to use for a matrix with distance-based disp. kernel
# Input:
# - D: a dispersal matrix with distance-based dispersal kernel
# - d: the desired mean dispersal rate
# Output: rescaled D matrix
rescale_D <- function(D, d=0.01) {
  D/(sum(D)/sum(D>0))*d
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





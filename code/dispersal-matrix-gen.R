# Generate random landscape with randomly drawn patch coordinates, each between 0 and 1.
# Input:
# - nPatch: Number of patches in the landscape.
# - dim: Landscape dimensionality (1 for a lakeshore; 2 for an area, etc.).
# Output:
# - An nPatch x dim data frame whose (i,j)th entry is the jth coordinate of patch i.
#   The rows are unnamed; the columns are named "V1", "V2", ..., "V<dim>".
randomCoords <- function(nPatch, dim = 2) {
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



#--------------------------------------------------------------------------------------
# Sample code:

# A 5-row, 2-column matrix of patch coordinates:
landscape <- randomCoords(5)
landscape

# 5x5 dispersal matrix with exponential kernel using those coordinates:
landscape |> dispMatrix(\(x) exp(-x))

# Three 5x5 dispersal matrices, with different dispersal kernels:
landscape |> dispMatrixList(list(\(x) exp(-x), \(x) 1/x, \(x) exp(-x^2)))

# The same, but organized into one large S*L x S*L (= 15 x 15) matrix:
landscape |> dispMatrixCommunity(list(\(x) exp(-x), \(x) 1/x, \(x) exp(-x^2)))

# Assuming exp. dispersal kernel and various dispersal distances, create dispersal matrix
# for 4 species (i.e., we get a 20 x 20 matrix):
landscape |> dispMatrixCommunityExp(c(3, 1, 0.5, 0.1))

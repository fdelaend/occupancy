library(deSolve)
library(rootSolve)
library(tidyverse)
extinctionThreshold <- 1e-3

# Matrix constructors -------------------------------------------------------
#' Build a block-diagonal matrix by repeating a base matrix.
#'
#' Creates a block-diagonal matrix with `p` copies of `A` placed on the diagonal.
#' Optionally introduces block-specific multiplicative variation in the entries
#' of `A`: each block is `A * V_i`, where `V_i` has i.i.d. Uniform(1 - vary, 1 + vary)
#' entries (element-wise scaling).
#'
#' @param A Numeric matrix. Template block placed along the diagonal.
#' @param p Integer >= 1. Number of blocks (i.e., number of copies of `A`).
#' @param vary Numeric >= 0. Magnitude of multiplicative variation applied
#'   independently to each entry of each block. `vary = 0` gives identical blocks.
#'   Note: values > 1 allow negative multipliers (since 1 - vary < 0).
#' @param ... Unused; included for interface compatibility.
#'
#' @return A numeric matrix of dimension `(nrow(A) * p) x (ncol(A) * p)` with
#'   `A`-derived blocks on the diagonal and zeros elsewhere.
#'
#' @details
#' - Variation is *multiplicative* and applied *element-wise* within each block.
#' - Off-diagonal blocks are exactly zero.
#' - Randomness comes from `runif()`. Set a seed externally for reproducibility.
#'
make_block_diagonal <- function(A, p, vary = 0, ...) {
  A_temp <- array(0, dim = dim(A) * p) # initialize with zeros
  for (i in 1:p) {
    # Block-specific multiplicative perturbation (element-wise)
    A_vary <- array(
      data = runif(prod(dim(A)), 1 - vary, 1 + vary),
      dim  = dim(A)
    )
    # Insert the perturbed block on the diagonal
    A_temp[(ncol(A) * i - ncol(A) + 1):(ncol(A) * i),
           (ncol(A) * i - ncol(A) + 1):(ncol(A) * i)] <- A * A_vary
  }
  return(A_temp)
}

#' Set the diagonal of a matrix.
#'
#' Replaces the diagonal elements of a matrix with specified values.
#' The input matrix is modified only on its diagonal; all off-diagonal
#' elements are left unchanged.
#'
#' @param A Numeric matrix. The matrix whose diagonal will be replaced.
#' @param d Numeric vector or scalar. Values assigned to the diagonal.
#'   If a scalar is provided, it is recycled to all diagonal elements.
#'
#' @return A numeric matrix with diagonal elements set to `d`.
#'
#' @details
#' - Length of `d` should be equal to `min(nrow(A), ncol(A))` or 1.
#' - No checks are performed on the length of `d`; recycling follows
#'   base R rules.
#'
set_diagonal <- function(A, d){
  diag(A) <- d
  A
}

#' Construct a species-specific dispersal matrix across locations.
#'
#' Builds a square dispersal matrix `D` describing exchange of individuals
#' among `p` locations for `n` species. The matrix is block-structured:
#' diagonal blocks (within-location dispersal) are zero, while off-diagonal
#' blocks represent dispersal between connected locations.
#'
#' For each connected pair of locations, dispersal occurs in both directions.
#' Off-diagonal blocks are diagonal matrices, allowing species-specific
#' dispersal rates but no cross-species dispersal.
#'
#' @param p Integer >= 1. Number of locations (patches).
#' @param n Integer >= 1. Number of species.
#' @param connectivity Numeric in [0, 1]. Fraction of all possible location
#'   pairs that are connected by dispersal.
#' @param d Numeric >= 0. Per-species dispersal rate between connected locations.
#'   The same rate is used for all species and all connections.
#' @param ... Unused; included for interface compatibility.
#'
#' @return A numeric matrix of dimension `(n * p) x (n * p)` representing
#'   species-specific dispersal among locations.
#'
#' @details
#' - Diagonal blocks are zero (no within-location dispersal).
#' - Off-diagonal blocks are diagonal matrices of dispersal rates.
#' - Connectivity is symmetric: if two locations are connected, dispersal
#'   occurs in both directions.
#' - Randomness arises from sampling which location pairs are connected.
#'   Set a random seed externally for reproducibility.
#'
make_D <- function(p = 3, n = 2, connectivity = 1, d = 1e-3, ...) {
  D <- array(0, dim = c(n * p, n * p))  # initialize dispersal matrix
  
  pos_links <- combn(c(1:p), 2)  # all possible location pairs
  nr_links  <- round(connectivity * ncol(pos_links), 0)  # number of connections
  links     <- pos_links[, sample(c(1:ncol(pos_links)), nr_links)]
  
  for (i in 1:ncol(links)) {
    # Species-specific dispersal matrices (no cross-species dispersal)
    D_12 <- diag(rep(d, n))  # dispersal from location 1 to 2
    D_21 <- diag(rep(d, n))  # dispersal from location 2 to 1
    
    pos_1 <- (n * links[1, i] - n + 1):(n * links[1, i])  # indices for location 1
    pos_2 <- (n * links[2, i] - n + 1):(n * links[2, i])  # indices for location 2
    
    D[pos_1, pos_2] <- D_12
    D[pos_2, pos_1] <- D_21
  }
  
  return(D)
}

#' Symmetrize a matrix while removing the diagonal.
#'
#' Constructs a symmetric matrix by combining the upper and lower triangular
#' parts of the input matrix. Diagonal elements are set to zero.
#'
#' Specifically, values from the upper triangle of `A` are copied to the lower
#' triangle (and vice versa), producing a symmetric matrix with no self-effects.
#'
#' @param A Numeric matrix. Input matrix to be symmetrized.
#'
#' @return A numeric symmetric matrix of the same dimension as `A`,
#'   with zeros on the diagonal.
#'
#' @details
#' - Only the strictly upper and lower triangular parts are used.
#' - Diagonal entries are discarded (set to zero).
#' - If `A` is already symmetric, this operation simply removes the diagonal.
#'
make_symmetric <- function(A) {
  A * upper.tri(A) + t(A) * lower.tri(A)
}

#' Construct a baseline interaction matrix with unit diagonal.
#'
#' Generates an `n x n` interaction matrix whose off-diagonal elements are drawn
#' from a normal distribution with mean `meanA` and standard deviation `sdA`.
#' Diagonal elements are set to 1, representing normalized self-regulation.
#'
#' @param meanA Numeric. Mean of the normal distribution used to generate
#'   off-diagonal interaction coefficients.
#' @param sdA Numeric >= 0. Standard deviation of the normal distribution.
#'   If `sdA = 0`, all off-diagonal elements are equal to `meanA`.
#' @param n Integer >= 1. Number of species.
#' @param ... Unused; included for interface compatibility.
#'
#' @return A numeric `n x n` matrix with diagonal entries equal to 1 and
#'   off-diagonal interaction coefficients drawn from `N(meanA, sdA^2)`.
#'
#' @details
#' - Diagonal entries represent intraspecific regulation and are fixed to 1.
#' - Off-diagonal entries represent interspecific interactions.
#' - Randomness arises from `rnorm()`. Set a random seed externally for
#'   reproducibility.
#'
make_A <- function(meanA = 0.6, sdA = 0, n = 3, ...) {
  matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>%
    set_diagonal(d = 1)
}

#' Generate spatially structured intrinsic growth rates for LV simulations.
#'
#' Constructs intrinsic growth-rate vectors for `n` species across `p` patches,
#' formatted for use in spatial Lotka–Volterra simulations. Growth rates are
#' generated such that species and patches are approximately equivalent on
#' average, while allowing controlled deviation from regional equivalence.
#'
#' @param n Integer >= 1. Number of species.
#' @param p Integer >= 1. Number of patches (locations).
#' @param k Numeric > 0. Strength of regional dominance. The species with the
#'   highest intrinsic growth rate is scaled by a factor `k` across all patches,
#'   while the remaining species are scaled by `1 / k`.
#' @param ... Unused; included for interface compatibility.
#'
#' @return A numeric vector of length `n * p` containing intrinsic growth rates,
#'   stacked by patch and species, suitable for direct use in spatial LV models.
#'
#' @details
#' - Initial growth-rate directions are sampled uniformly from the positive
#'   orthant of the unit hypersphere.
#' - An iterative rescaling procedure enforces approximate equality of mean
#'   growth rates across species and across patches.
#' - Regional equivalence is then relaxed by amplifying one species’ growth
#'   rate across all patches via the parameter `k`.
#' - Randomness arises from `rnorm()`. Set a random seed externally for
#'   reproducibility.
#'
make_R_spatial <- function(n, p, k = 1, ...) {
  # Sample directions from the positive orthant of a hypersphere
  rawRs <- abs(matrix(rnorm(p * n), nrow = p))
  rawRs <- t(rawRs / sqrt(rowSums(rawRs^2)))
  
  # Enforce approximate equality of mean growth rates across
  # species (columns) and patches (rows) via iterative scaling
  for (iteration in 1:2) {
    rawRs <- rawRs %*% diag(1 / colMeans(rawRs))
    rawRs <- diag(1 / rowMeans(rawRs)) %*% rawRs
  }
  
  # Relax regional equivalence by making one species globally dominant
  rawRs <- t(rawRs) %*% diag(c(k, rep(1 / k, n - 1)))
  colnames(rawRs) <- 1:n
  
  # Return growth rates as a stacked vector (patch × species)
  Rs <- as_tibble(rawRs) %>%
    mutate(location = 1:p) %>%
    pivot_longer(!location)
  
  return(Rs$value)
}

#' Generate random spatial coordinates for habitat patches.
#'
#' Creates a random landscape by assigning each of `nPatch` habitat patches
#' a position in a `dim`-dimensional unit space. Patch coordinates are drawn
#' independently from a uniform distribution on [0, 1].
#'
#' @param nPatch Integer >= 1. Number of habitat patches.
#' @param dim Integer >= 1. Spatial dimensionality of the landscape
#'   (e.g., 1 for a linear landscape, 2 for a planar landscape).
#'
#' @return A data frame with `nPatch` rows and `dim` columns, where each row
#'   gives the spatial coordinates of a patch. Columns are named `V1`, `V2`,
#'   ..., `V<dim>`.
#'
#' @details
#' - Coordinates are drawn independently for each dimension.
#' - The spatial domain is the unit hypercube [0, 1]^dim.
#' - Randomness arises from `runif()`. Set a random seed externally for
#'   reproducibility.
#'
make_randomCoords <- function(nPatch, dim = 2) {
  matrix(runif(nPatch * dim, 0, 1), nrow = nPatch, ncol = dim) %>%
    as.data.frame(row.names = NULL)
}

## Dispersal kernel functions --------------------------------------------------------
#' Construct a distance-based dispersal matrix from patch coordinates.
#'
#' Computes pairwise distances between habitat patches and converts them into
#' dispersal rates using a user-specified dispersal kernel. The resulting matrix
#' represents directed dispersal among patches, with self-dispersal excluded.
#'
#' @param coords Numeric matrix, data frame, or tibble. Patch coordinates, with
#'   patches in rows and spatial dimensions in columns.
#' @param kernel Function. Dispersal kernel mapping distance to dispersal rate.
#'   Must accept a single numeric argument (distance) and return a numeric value
#'   or vector of the same length.
#'
#' @return A numeric matrix where entry `(i, j)` gives the dispersal rate from
#'   patch `j` to patch `i`. Diagonal entries are set to zero.
#'
#' @details
#' - Pairwise Euclidean distances are computed using `dist()`.
#' - The kernel is applied element-wise to the distance matrix.
#' - The resulting matrix is symmetric if the kernel is symmetric.
#' - Self-dispersal is excluded by zeroing the diagonal.
#'
dispMatrix <- function(coords, kernel) {
  coords %>%                     # patch coordinates
    dist() %>%                   # pairwise Euclidean distances
    as.matrix() %>%              # convert to full distance matrix
    kernel() %>%                 # apply dispersal kernel
    unname() %>%                 # remove row and column names
    (`diag<-`)(0)                # remove self-dispersal
}

#' Construct species-specific dispersal matrices for a shared landscape.
#'
#' Generates a list of dispersal matrices from a common set of patch coordinates,
#' allowing each species to have its own distance-dependent dispersal kernel.
#' All matrices share the same spatial structure but differ in the functional
#' form of dispersal.
#'
#' @param coords Matrix, data frame, or tibble. Patch coordinates with patches
#'   in rows and spatial dimensions in columns.
#' @param kernelList List of functions. Each function specifies a dispersal
#'   kernel for a species and must accept a numeric vector or matrix of distances
#'   and return values of the same shape.
#'
#' @return A list of numeric square matrices of dimension
#'   `nrow(coords) x nrow(coords)`. The `i`th list element corresponds to the
#'   dispersal matrix of species `i`.
#'
#' @details
#' - All species disperse over the same landscape geometry.
#' - Species differ only in the functional form of their dispersal kernels.
#' - Diagonal entries of each matrix are zero (no self-dispersal), as enforced
#'   by `dispMatrix()`.
#'
dispMatrixList <- function(coords, kernelList) {
  # Apply the dispersal-matrix constructor to each species-specific kernel
  Map(\(f) dispMatrix(coords, f), kernelList)
}

#' Creates a square matrix whose entries are all zero except for a single
#' diagonal element, which is set to one. 
#' 
#' This matrix is useful for constructing
#' block-structured or community-level matrices via Kronecker products.
#'
#' @param dim Integer >= 1. Dimension of the square matrix.
#' @param n Integer. Index of the diagonal entry to be set to one.
#'   Must satisfy `1 <= n <= dim`.
#'
#' @return A numeric `dim x dim` matrix with a single one at position `(n, n)`
#'   and zeros elsewhere.
#'
#' @details
#' - The resulting matrix acts as a selector for the `n`th diagonal component.
#' - Commonly used in Kronecker-sum constructions of large block matrices.
#'
oneHotMatrix <- function(dim, n) {
  rep(0, dim) |>          # initialize zero vector
    replace(n, 1) |>      # set the nth entry to one
    diag()                # place vector on the diagonal
}

#' Construct the full community-wide dispersal matrix.
#'
#' Builds an expanded dispersal matrix of dimension `(S * L) x (S * L)`,
#' where `S` is the number of species and `L` the number of patches.
#' Each species disperses over the same landscape but may have its own
#' distance-dependent dispersal kernel.
#'
#' The resulting matrix is assembled using Kronecker products, embedding each
#' species-specific `L x L` dispersal matrix into the appropriate block of
#' the full community matrix.
#'
#' @param coords Matrix, data frame, or tibble. Patch coordinates with patches
#'   in rows and spatial dimensions in columns.
#' @param kernelList List of functions. Each function specifies a dispersal
#'   kernel for one species and must accept a numeric vector or matrix of
#'   distances and return values of the same shape.
#'
#' @return A numeric square matrix of dimension `(S * L) x (S * L)` representing
#'   community-wide dispersal across species and patches.
#'
#' @details
#' - Let `D_i` denote the `L x L` dispersal matrix of species `i`.
#' - The full matrix is constructed as:
#'     sum_i ( D_i ⊗ E_i ),
#'   where `E_i` is a one-hot `S x S` diagonal matrix selecting species `i`,
#'   and `⊗` denotes the Kronecker product.
#' - Species disperse independently (no cross-species dispersal).
#'
dispMatrixCommunity <- function(coords, kernelList) {
  S <- length(kernelList)         # number of species
  L <- nrow(coords)               # number of patches
  
  dispMatrices <- coords %>% 
    dispMatrixList(kernelList)    # list of S dispersal matrices (each L x L)
  
  fullDisp <- matrix(0, S * L, S * L)  # initialize community matrix
  
  # Embed each species' dispersal matrix using Kronecker products
  for (i in 1:S)
    fullDisp <- fullDisp + dispMatrices[[i]] %x% oneHotMatrix(S, i)
  
  fullDisp
}

#' Construct a community-wide dispersal matrix with exponential kernels.
#'
#' Convenience wrapper around `dispMatrixCommunity()` that assumes all species
#' disperse according to exponential distance kernels. Instead of supplying a
#' list of kernel functions, the user provides a vector of characteristic
#' dispersal distances, one per species.
#'
#' @param coords Matrix, data frame, or tibble. Patch coordinates with patches
#'   in rows and spatial dimensions in columns.
#' @param dispDistanceVector Numeric vector of length `S`, where `S` is the
#'   number of species. The `i`th entry specifies the characteristic dispersal
#'   distance of species `i`.
#'
#' @return A numeric square matrix of dimension `(S * L) x (S * L)`, where
#'   `L = nrow(coords)`, representing community-wide dispersal with
#'   species-specific exponential kernels.
#'
#' @details
#' - For species `i`, the dispersal kernel is:
#'     \( K_i(x) = \exp(-x / d_i) \),
#'   where `d_i` is the corresponding entry of `dispDistanceVector`.
#' - Species differ only in their characteristic dispersal distance.
#' - The full community matrix is assembled using the same Kronecker-product
#'   construction as in `dispMatrixCommunity()`.
#'
dispMatrixCommunityExp <- function(coords, dispDistanceVector) {
  # Build a list of species-specific exponential dispersal kernels
  Map(\(d) \(x) exp(-x / d), dispDistanceVector) |>
    # Construct the community-wide dispersal matrix
    dispMatrixCommunity(coords = coords)
}

#' Rescale a dispersal matrix to achieve a target mean dispersal rate.
#'
#' Rescales a dispersal matrix derived from a distance-based dispersal kernel
#' so that the mean of its positive (non-zero) entries equals a specified
#' target value.
#'
#' @param D Numeric matrix. Dispersal matrix whose non-zero entries represent
#'   distance-dependent dispersal rates.
#' @param d Numeric > 0. Desired mean dispersal rate across all non-zero
#'   entries of `D`.
#'
#' @return A numeric matrix of the same dimension as `D`, rescaled so that
#'   `mean(D[D > 0]) = d`.
#'
#' @details
#' - Only non-zero entries are used to compute the mean dispersal rate.
#' - Zero entries (e.g., self-dispersal or unconnected patch pairs) remain zero.
#' - This operation preserves the relative structure imposed by the dispersal
#'   kernel and rescales only its overall magnitude.
#'
rescale_D <- function(D, d = 0.01) {
  D / (sum(D) / sum(D > 0)) * d
}

## LV dynamics -------------------------------------

#' Smooth cutoff function for state variables.
#'
#' Applies a smooth, bounded cutoff to state variables to improve numerical
#' stability when solving ordinary differential equations. 
#'
#' @param n Numeric vector. State variable values to be constrained.
#'
#' @return Numeric vector of the same length as `n`, with values smoothly
#'   mapped to the interval [0, 1].
#'
#' @details
#' - For `0 < n < 1`, values are transformed using a cubic polynomial chosen
#'   to ensure smoothness at the boundaries.
#' - Values `n >= 1` are set to 1.
#' - Values `n <= 0` are mapped to 0.
#' - The function is continuous and has continuous first derivatives,
#'   avoiding numerical instabilities associated with hard thresholds.
#'
cutoff <- function(n) {
  ifelse(
    n < 1,
    (1 * (n > 0)) * (n * n * n * (10 + n * (-15 + 6 * n))),
    1
  )
}

#' Lotka–Volterra population dynamics with numerical stabilization.
#'
#' Implements a continuous-time Lotka–Volterra (LV) model for interacting
#' species.
#' Intrinsic growth rates and interaction coefficients are supplied via
#' the parameter list, while a smooth cutoff is applied to enhance numerical
#' stability at very low abundances.
#'
#' @param time Numeric. Time variable (required by ODE solvers; not used
#'   explicitly in the model).
#' @param state Numeric vector. Current state vector of species abundances.
#' @param pars List of model parameters with components:
#'   \describe{
#'     \item{R}{Numeric vector of intrinsic growth rates.}
#'     \item{A}{Numeric interaction matrix defining linear species interactions.}
#'   }
#'
#' @return A list containing the time derivatives of `state`, formatted for
#'   use with ODE solvers such as `deSolve::ode()`.
#'
#' @details
#' - The model has the form:
#'   \( \dot{N} = N \odot (R + A N) \),
#'   where `⊙` denotes element-wise multiplication.
#' - A smooth cutoff function is applied to suppress numerical artefacts when
#'   abundances approach zero.
#' - The cutoff affects only extremely small values and does not alter the
#'   qualitative dynamics of the LV system.
#'
run_LV <- function(time, state, pars) {
  list(as.numeric(
    (state * (pars$R + pars$A %*% state)) * cutoff(state / 1e-10)
  ))
}

#' Spatial Lotka–Volterra dynamics with dispersal.
#'
#' Implements a spatially explicit Lotka–Volterra (LV) model in which species
#' interact locally and disperse among patches. Local population dynamics are
#' governed by intrinsic growth rates and linear interaction terms, while
#' dispersal is represented by a community-wide dispersal matrix.
#'
#' @param time Numeric. Time variable (required by ODE solvers; not used
#'   explicitly in the model).
#' @param state Numeric vector. Current state vector of species abundances,
#'   stacked across patches and species.
#' @param pars List of model parameters with components:
#'   \describe{
#'     \item{R}{Numeric vector of intrinsic growth rates for all species
#'       in all patches.}
#'     \item{A}{Numeric interaction matrix defining linear species interactions.}
#'     \item{D}{Numeric dispersal matrix describing exchange among patches
#'       and species.}
#'   }
#'
#' @return A list containing the time derivatives of `state`, formatted for
#'   use with ODE solvers such as `deSolve::ode()`.
#'
#' @details
#' - The local (within-patch) dynamics follow:
#'   \( \dot{N} = N \odot (R - A N) \).
#' - Dispersal contributes two terms: emigration (`- N \odot \sum_j D_{j\cdot}`)
#'   and immigration (`D N`).
#' - The subtraction of `colSums(D)` from `R` ensures mass balance by accounting
#'   for total emigration rates.
#' - The model assumes no cross-species dispersal unless encoded in `D`.
#'
run_LV_spatial <- function(time, state, pars) {
  list(
    (state * ((pars$R - colSums(pars$D)) - drop(pars$A %*% state))) +
      drop(pars$D %*% state)
  )
}

#' Jacobian of the spatial Lotka–Volterra model.
#'
#' Computes the Jacobian matrix of the spatially explicit Lotka–Volterra (LV)
#' system with dispersal, evaluated at the current state. This matrix is used
#' for local stability analysis and by ODE solvers that exploit analytic
#' derivatives.
#'
#' @param time Numeric. Time variable (required by ODE solvers; not used
#'   explicitly in the Jacobian).
#' @param state Numeric vector. Current state vector of species abundances,
#'   stacked across patches and species.
#' @param pars List of model parameters with components:
#'   \describe{
#'     \item{R}{Numeric vector of intrinsic growth rates.}
#'     \item{A}{Numeric interaction matrix defining linear species interactions.}
#'     \item{D}{Numeric dispersal matrix describing exchange among patches
#'       and species.}
#'   }
#'
#' @return A numeric square matrix giving the Jacobian of the spatial LV system
#'   evaluated at `state`.
#'
#' @details
#' - The Jacobian corresponds to the linearization of:
#'   \( \dot{N} = N \odot (R - A N) + D N - N \odot \sum_j D_{j\cdot} \).
#' - The first term represents density-dependent local growth.
#' - The second term (`- state * A`) captures how interactions affect growth
#'   rates through changes in species abundances.
#' - The dispersal matrix `D` contributes additively to the Jacobian.
#'
Jac_LV_spatial <- function(time, state, pars) {
  (pars$R - colSums(pars$D) - drop(pars$A %*% state)) * diag(length(pars$R)) -
    state * pars$A + pars$D
}

#' Compute equilibrium (steady state) of the spatial LV model.
#'
#' Integrates the spatial Lotka–Volterra system with dispersal until it reaches
#' a steady state, using `deSolve::runsteady()`. This is used to obtain the
#' equilibrium community composition from an initial condition.
#'
#' @param n Integer >= 1. Number of species.
#' @param R Numeric vector. Intrinsic growth rates for all species across all
#'   patches (length `n * p`).
#' @param A Numeric matrix. Interaction matrix (as in `run_LV_spatial()`).
#' @param D Numeric matrix. Dispersal matrix (as in `run_LV_spatial()`).
#' @param max_time Numeric. (Currently unused) Maximum simulation time. Included
#'   for interface compatibility with alternative ODE-based implementations.
#' @param N0 Numeric vector. Initial abundances (length `n * p`).
#' @param tstep Numeric. (Currently unused) Time step used only in the commented
#'   `ode()` alternative.
#' @param ... Unused; included for interface compatibility.
#'
#' @return A numeric vector `y` giving the steady-state abundances (length `n * p`).
#'   Note: this is the `$y` element returned by `deSolve::runsteady()`.
#'
#' @details
#' - The number of patches `p` is inferred as `length(R) / n`.
#' - State names are set to `"species_patch"` for readability.
#' - The Jacobian `Jac_LV_spatial()` is supplied to improve solver performance.
#'
get_dynamics <- function(n, R, A, D, max_time = 200, N0, tstep = 10, ...) {
  p <- length(R) / n  # infer number of patches
  
  # Name state variables as "<species>_<patch>" for readability
  names(N0) <- paste(rep(1:n, p), "_", rep(1:p, each = n), sep = "")
  
  # Integrate until equilibrium
  runsteady(
    y       = N0,
    func    = run_LV_spatial,
    parms   = list(R = R, A = A, D = D),
    jacfunc = Jac_LV_spatial,
    mf      = 21,
    stol    = 1e-12
  )$y
  
}

#' Compute and return equilibrium abundances of the spatial LV model.
#'
#' Convenience wrapper around `get_dynamics()` that extracts and formats
#' the steady-state solution of the spatial Lotka–Volterra system.
#'
#' @param n Integer >= 1. Number of species.
#' @param R Numeric vector. Intrinsic growth rates (length `n * p`).
#' @param A Numeric matrix. Interaction matrix (as in `run_LV_spatial()`).
#' @param D Numeric matrix. Dispersal matrix (as in `run_LV_spatial()`).
#' @param max_time Numeric. Passed to `get_dynamics()` (currently unused in
#'   steady-state computation).
#' @param N0 Numeric vector. Initial abundances (length `n * p`).
#' @param ... Unused; included for interface compatibility.
#'
#' @return Formatted equilibrium abundances as produced by `organise_ode()`.
#'
#' @details
#' - Calls `get_dynamics()` to compute the steady state using `runsteady()`.
#' - Applies `organise_ode()` to restructure the output into a tidy format.
#' - This function serves as the primary entry point for obtaining equilibrium
#'   abundances in spatial simulations.
#'
get_NHat <- function(n, R, A, D, max_time = 200, N0, ...) {
  get_dynamics(n = n, R = R, A = A, D = D,
               max_time = max_time, N0 = N0) %>%
    organise_ode()
}

#' Reorganize ODE output into tidy species–patch format.
#'
#' Reshapes the output of a spatial Lotka–Volterra simulation into a tidy
#' data structure with explicit species and location identifiers. This
#' function assumes that state variables are named using the pattern
#' `"species_patch"`.
#'
#' @param dynamics Numeric vector or object containing equilibrium abundances,
#'   with names encoding species and patch indices as `"sp_location"`.
#' @param ... Unused; included for interface compatibility.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{sp}{Integer species index.}
#'     \item{location}{Integer patch index.}
#'     \item{value}{Equilibrium abundance (density).}
#'   }
#'
#' @details
#' - Uses `enframe()` to convert named values into a tidy format.
#' - Species and patch indices are extracted from the state names.
#' - This function is typically applied to the output of `get_dynamics()`
#'   or `get_NHat()`.
#'
organise_ode <- function(dynamics, ...) {
  enframe(dynamics, name = "sp_location") %>%
    separate_wider_delim(
      sp_location,
      delim = "_",
      names = c("sp", "location")
    ) %>%
    mutate(
      sp = as.integer(sp),
      location = as.integer(location)
    )
}

#' Summarize equilibrium outcomes of the spatial LV model.
#'
#' Computes patch-level and richness-level summaries from equilibrium
#' abundances of the spatial Lotka–Volterra model. Species are classified
#' as present or extinct based on a fixed extinction threshold.
#'
#' @param NHat Tibble. Equilibrium abundances as returned by `organise_ode()`,
#'   with columns `sp`, `location`, and `value`.
#' @param R Numeric vector. Intrinsic growth rates corresponding to each
#'   species–patch combination.
#' @param n Integer >= 1. Total number of species in the regional pool.
#' @param ... Unused; included for interface compatibility.
#'
#' @return A tibble summarizing equilibrium properties by local species
#'   richness, with columns:
#'   \describe{
#'     \item{m}{Local species richness (number of persisting species).}
#'     \item{nrPatches}{Number of patches exhibiting richness `m`.}
#'     \item{NTotal}{Mean total biomass per patch at richness `m`.}
#'     \item{meanRPer}{Mean intrinsic growth rate of persisting species,
#'       averaged across patches with richness `m`.}
#'   }
#'
#' @details
#' - Species are considered present if their equilibrium abundance exceeds
#'   `extinctionThreshold`.
#' - Local richness is computed per patch as the number of persisting species.
#' - Results are first summarized at the patch level and then aggregated
#'   across patches with identical richness.
#' - Richness levels with zero patches are retained for completeness.
#'
summarize_ode <- function(NHat, R, n, ...) {
  NHat %>%
    mutate(
      present = value > extinctionThreshold,
      R = R
    ) %>%
    summarize(
      m         = as_factor(sum(present)),
      NTotal    = sum(value),
      meanRPer  = sum(R * present) / sum(present),
      .by = location
    ) %>%
    mutate(
      m = fct_expand(m, as.character(c(1:n)))
    ) %>%
    group_by(m, .drop = FALSE) %>%
    summarize(
      nrPatches = n(),            # number of patches with m species
      NTotal    = mean(NTotal),   # mean total biomass per patch
      meanRPer  = mean(meanRPer)  # mean intrinsic growth rate of persisting species
    )
}



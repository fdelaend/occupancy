
source("tools-simulations.R")
parallel_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# SIMULATIONS ----
## Simulations ------
Sims <- expand_grid(n = c(6), meanA = c(0.2, 0.4, 0.8, 1), #, 6; 0.4, 0.4, 0.6, 
                    d = seq(-6,-4, length.out=6), vary=c(0, 0.1), k=c(1, 1.5),
                    cvA = c(0, 0.2), p = c(10, 50, 100), rep = parallel_id) |> #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters: d, sdA
  mutate(d = 10^d,
         sdA = cvA * meanA) |>
  #Construct symmetric A-matrix
  (\(x) mutate(x, A = pmap(x, \(meanA, sdA, n, p, vary, ...) 
                           matrix(data = rnorm(n^2, meanA, sdA), ncol = n) |>
                             make_symmetric() |> 
                             set_diagonal(d=1) |>
                             make_block_diagonal(p=p, var=vary))))() |>
  #Sample growth rates
  (\(x) mutate(x, R = pmap(x, make_R_spatial)))() |>
  #Make regular D matrix (representing well-mixed spatial network)
  (\(x) mutate(x, regularD = pmap(x, make_D)))() |>
  #Make spatially-explicity D matrix (with exponential dispersal kernel)
  mutate(coords = map(p, ~make_randomCoords(nPatch=.x)), #generate coordinates for every patch
         distances = map(n, ~seq(4,0.5, length.out=.x))) |> #generate characteristic distances of the species
  mutate(exponentialD = map2(coords, distances, ~dispMatrixCommunityExp(coords=.x, #Make a D matrix for exponential decay
                                                                        dispDistanceVector=.y)),
         exponentialD = map2(exponentialD, d, ~rescale_D(D=.x, d=.y))) |> #and rescale to have same mean D as in regular case
  mutate(N0 = map(R, ~.x/10)) |> #set initial conditions equal to carrying capacity w/o interactions or dispersal, divided by 10
  pivot_longer(cols = c("regularD","exponentialD"), names_to = "dispType", values_to = "D") |> # pivot the two sorts of D matrices to make them a factor
  #Simulate the network
  (\(x) mutate(x, NHat = pmap(x, get_NHat)))() |>
  #1/ Summarize the simulated data: per m, compute the nr of patches and total biomass of an average patch
  (\(x) mutate(x, summaryM = pmap(x, \(NHat, R, n,...) 
                                  NHat |>
                                    mutate(present = density>extinctionThreshold,
                                           R=R) |>
                                    summarize(m = as_factor(sum(present)),
                                              NTotal = sum(density),
                                              meanRPer = sum(R*present)/sum(present), 
                                              .by = location) |>
                                    mutate(m = fct_expand(m, as.character(c(1:n)))) |>
                                    group_by(m, .drop=F) |>
                                    summarize(nrPatches = n(),#nr of patches with m sp.
                                              NTotal = mean(NTotal),#total biomass in a patch with m sp.
                                              meanRPer = mean(meanRPer)))))()  |>#mean r of persisting sp.
  #2/proportion of patches in which all n species persist
  mutate(propPatchesN = 1/p*map2_dbl(summaryM, n, ~ (.x |> filter(m==.y))$nrPatches)) |>
  #3/total density across all patches of a species
  mutate(NTotalK = map(NHat, ~ .x |> 
                         summarize(NTotalK = sum(density), .by = sp)))

saveRDS(Sims, file=paste("../simulated-data/",parallel_id,"data.RDS",sep=""))

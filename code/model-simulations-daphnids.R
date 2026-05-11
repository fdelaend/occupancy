
source("tools-simulations.R")
parallel_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') # Get task ID from the cluster

# SIMULATIONS ----
## Simulations ------
Sims <-
  expand_grid(n = c(3), meanA = c(1, 1.2, 1.4), #preemptive competition, representative for Daphnid system
              d = c(seq(-6, -4, length.out=6)), #dispersal rate estimations based on Dubart
              vary=c(0, 0.1, 0.2, 0.3), k=c(1, 1.2, 1.4), #uncertain about spatial variation of interactions or regional inequality
              cvA = c(1e-2, 0.2), p = seq(3, 100, 5), # uncertain about the CV of species interactions; nr of patches as in Daphnid system
              rep = parallel_id) %>% 
  #Make parameters: d, sdA
  mutate(d = 10^d,
         sdA = cvA * meanA) %>%
  #Construct symmetric A-matrix
  (\(x) mutate(x, A = pmap(x, \(meanA, sdA, n, p, vary, ...) {
                           matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>%
                             make_symmetric() %>%
                             set_diagonal(d=1) %>%
                             make_block_diagonal(p=p, var=vary)
  } )))() %>%
  #Sample growth rates
  (\(x) mutate(x, R = pmap(x, make_R_spatial)))() %>%
  #Make regular D matrix (representing well-mixed landscape)
  (\(x) mutate(x, regularD = pmap(x, make_D)))() %>%
  #Make spatially-explicity D matrix (with exponential dispersal kernel)
  mutate(coords = map(p, ~make_randomCoords(nPatch=.x)), #generate coordinates for every patch, representative for Daphnid system
         distances = map(n, ~seq(0.1,0.5, length.out=.x))) %>% #generate characteristic distances of the species; representative for clustered dispersal
  mutate(exponentialD = map2(coords, distances, ~dispMatrixCommunityExp(coords=.x, #Make a D matrix for exponential decay
                                                                        dispDistanceVector=.y)),
         exponentialD = map2(exponentialD, d, ~rescale_D(D=.x, d=.y))) %>% #and rescale to have same mean D as in regular case
  mutate(N0 = map(R, ~.x/10)) %>% #set initial conditions equal to carrying capacity w/o interactions or dispersal, divided by 10
  pivot_longer(cols = c("regularD","exponentialD"), names_to = "dispType", values_to = "D") %>% # pivot the two sorts of D matrices to make them a factor
  #Simulate the landscape
  (\(x) mutate(x, NHat = pmap(x, get_NHat, .progress = TRUE)))() #, .progress = TRUE

## Summarize results ----

results <- Sims %>%
  #1/ Summarize the simulated data: per m, compute the nr of patches and total biomass of an average patch
  (\(x) mutate(x, summaryM = pmap(x, summarize_ode)))() %>%
  #2/proportion of patches in which all n species persist
  mutate(propPatchesN = 1/p*map2_dbl(summaryM, n, ~ (.x %>% filter(m==.y))$nrPatches),
         propPatches2 = 1/p*map_dbl(summaryM, ~ (.x |> filter(m==2))$nrPatches)) %>%
  #3/total density across all patches of a species
  mutate(NTotalK = map(NHat, ~ .x %>%
                         summarize(NTotalK = sum(value), 
                                   .by = sp))) %>%
  select(!A & !D & !distances & !N0 & !coords) #ditch all unneeded data

## Save results ----
write_rds(results, 
          file=str_c("../simulated-data/", 
                     parallel_id,"data.RDS"))

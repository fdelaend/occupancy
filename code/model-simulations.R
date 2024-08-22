
source("tools-simulations.R")
parallel_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# SIMULATIONS ----
## Simulations ------
Sims <- expand_grid(n = c(6), meanA = c(0.2, 0.4, 0.8, 1), #, 6; 0.4, 0.4, 0.6, 
                    d = seq(-6,-4, length.out=6), vary=c(0, 0.1), k=c(1, 1.5),
                    cvA = c(0, 0.2), p = c(10, 100), rep = parallel_id) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters
  mutate(d = 10^d,
         sdA = cvA * meanA) %>%
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% #make symmetric and ditch diagonal; set self to 1
      make_symmetric() %>% set_diagonal(d=1))) %>% 
  mutate(A = pmap(., make_block_diagonal)) %>% #make A spatial
  mutate(R = pmap(., make_R_spatial)) %>% #make spatial R (note that the emigration is subtracted later so this is the real local R)
  mutate(regularD = pmap(., make_D)) %>% #make classic dispersal matrix
  mutate(coords = map(p, ~make_randomCoords(nPatch=.x)), #generate coordinates for every patch
         distances = map(n, ~seq(4,0.5, length.out=.x))) %>% #generate characteristic distances of the species
  mutate(exponentialD = map2(coords, distances, ~dispMatrixCommunityExp(coords=.x, #Make a D matrix for exponential decay
                                                                        dispDistanceVector=.y)),
         exponentialD = map2(exponentialD, d, ~rescale_D(D=.x, d=.y))) %>%
  mutate(N0 = map(R, ~.x/10)) %>% #set initial conditions equal to carrying capacity w/o interactions or dispersal, divided by 10
  pivot_longer(cols = c("regularD","exponentialD"), names_to = "dispType", values_to = "D") %>% # pivot the two sorts of D matrices to make them a factor
  #Simulate the network
  mutate(NHat = pmap(., get_NHat)) %>%
  #1/ Summarize: per m, compute the nr of patches and total biomass of an average patch
  mutate(summaryM = pmap(., function(NHat, R, n,...) {
    NHat %>% mutate(present = density>extinctionThreshold) %>%
      mutate(R=R) %>%
      group_by(location) %>%
      summarize(m = as_factor(sum(present)),
                NTotal = sum(density),
                meanRPer = sum(R*present)/sum(present)) %>%
      mutate(m = fct_expand(m, as.character(c(1:n)))) %>%
      group_by(m, .drop=F) %>%
      summarize(nrPatches = n(),#nr of patches with m sp.
                NTotal = mean(NTotal),#total biomass in a patch with m sp.
                meanRPer = mean(meanRPer)) %>% #mean r of persisting sp.
      ungroup()})) %>%
  #2/proportion of patches in which all n species persist
  mutate(propPatchesN = 1/p*map2_dbl(summaryM, n, ~ (.x %>% filter(m==.y))$nrPatches)) %>%
  #3/total density across all patches of a species
  mutate(NTotalK = pmap(., function(NHat,...) {
    NHat %>% 
      group_by(sp) %>% #sum across patches for every sp.
      summarize(NTotalK = sum(density))})) #so you can check it's comparable across sp

saveRDS(Sims, file=paste("../data/",parallel_id,"data.RDS",sep=""))

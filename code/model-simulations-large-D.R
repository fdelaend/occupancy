
# Illustration of what happens when dispersal is large

Sims <- expand_grid(n = c(6), meanA = c(0.8), #, 6; 0.4, 0.4, 0.6, 
                    d = seq(-3,0, length.out=5), vary=0, k=c(1, 1.5),
                    cvA = 0.2, p = c(10, 20, 30, 40, 50), rep = 1) |> #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
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
  (\(x) mutate(x, D = pmap(x, make_D)))() |>
  mutate(N0 = map(R, ~.x/10)) |> #set initial conditions equal to carrying capacity w/o interactions or dispersal, divided by 10
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
                         summarize(NTotalK = sum(density), .by = sp))) |>
  select(!A & !N0) #ditch all unneeded data

ggplot(Sims) + 
  aes(x=p, y=log10(d), col=propPatchesN)+
  geom_point() +
  facet_grid(k~.)
  
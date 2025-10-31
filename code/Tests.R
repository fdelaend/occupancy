
Sims <- expand_grid(n = c(200), meanA = c(0.4), it = c(1:2), #, 6; 0.4, 0.4, 0.6, 
                    d = c(-20), vary=c(0), k=c(1),
                    sdA = c(1), p = c(3), rep = parallel_id) |> #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  mutate(d = 10^d, 
         meanA = meanA / n, #scaling
         sdA = sdA/sqrt(n)) |> #scaling
  #Construct A-matrix
  (\(x) mutate(x, A = pmap(x, \(meanA, sdA, n, p, vary, ...) 
                           matrix(data = rnorm(n^2, meanA, sdA), ncol = n) |>
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
  filter(dispType == "regularD") |>
  (\(x) mutate(x, dynamics = pmap(x, get_dynamics)))() |>
  mutate(limit_index = map_dbl(dynamics, ~get_limit_index(.x, tail_size = 100))) |>
  (\(x) mutate(x, dynamics = pmap(x, organise_ode)))()
  
#Plot the dynamics
Sims_plot <- Sims |>
  select(n, meanA, d, vary, k, sdA, p, it,
         dynamics, limit_index) |>
  unnest(dynamics) |>
  filter(vary == 0, time %in% c(1:200))

ggplot(Sims_plot) + 
  theme_bw() +
  geom_line(aes(x = time, y = N, colour = as.factor(sp),
                group = interaction(location, sp)),
            show.legend = F) +
  geom_text(aes(x = 50, y = 3, label = limit_index)) +
  facet_grid(it ~ sdA)



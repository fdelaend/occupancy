
source("tools-simulations.R")

Sims <-
  expand_grid(n = c(10), meanA = c(0.2, 0.5), #, 6; 0.4, 0.4, 0.6,
              vary=c(0), k=c(1), cvA = c(1e-2), 
              p = 50, rep = c(1)) |> #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters: sdA
  mutate(sdA = cvA * meanA) |>
  #Construct symmetric A-matrix
  (\(x) mutate(x, A = pmap(x, \(meanA, sdA, n, p, vary, ...) {
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) |>
      make_symmetric() |>
      set_diagonal(d=1) |>
      make_block_diagonal(p=p, var=vary)
  } )))() |>
  #Sample growth rates
  (\(x) mutate(x, R = pmap(x, make_R_spatial)))() |>
  expand_grid(d = seq(0, 0.075, length.out=30)) |>
  #Make regular D matrix (representing well-mixed spatial network)
  (\(x) mutate(x, D = pmap(x, make_D)))() |>
  mutate(N0 = map(R, ~.x/10)) |> #set initial conditions equal to carrying capacity w/o interactions or dispersal, divided by 10
  #Simulate the network
  (\(x) mutate(x, NHat = pmap(x, get_NHat, .progress = TRUE)))() #, .progress = TRUE

Sims_plot <- Sims |>
  select(n, meanA, p, rep, d, NHat) |>
  unnest(NHat) |>
  summarise(mean = mean(value),
            sd = sd(value),
            .by = c(sp, d, meanA))

ggplot(Sims_plot) + 
  aes(x = d, y = mean, col = factor(sp)) + 
  theme_bw() +
  #geom_errorbar(width = 0.0001) + 
  geom_point() +
  facet_grid(meanA~., 
             labeller = label_bquote(cols=paste("a = ",.(meanA)))) +
  ylim(0, 0.6) + 
  labs(y = "mean density across patches", 
       col = "species") + 
  geom_vline(xintercept = 0.005, lty = "dashed")

ggsave(paste0("../figures/d.pdf"), width=5, height = 5, device = "pdf")  
  
  

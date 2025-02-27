
test <- tibble(N0 = rnorm(1000, 0.1, sd=1)) |>
  mutate(N0=if_else(N0>0, N0, 0)) |>
  expand_grid(A = 10^(seq(0.6, 1.4, 0.1)), 
              d = 2*10^-2) |>
  mutate(N = N0 + d*A) |>
  summarise(richness = sum(N>1), .by = A)
  
ggplot(test) + 
  aes(x = N, col = as.factor(A)) + 
  geom_density() + 
  geom_vline(xintercept = 1)

ggplot(test) + 
  aes(x = log10(A), y = log10(richness)) + 
  geom_point() 

lm(log10(richness) ~ log10(A), data = test)

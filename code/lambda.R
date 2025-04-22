
expand.grid(p = c(5:50)) |>
  mutate(lambda = map_dbl(p, ~max(abs(eigen(make_D(.x))$values)))) |>
  ggplot() + 
  aes(x = p, y = lambda) + 
  geom_point()
  
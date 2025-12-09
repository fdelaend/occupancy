
test <- counts_des |>
  filter(richness > 0, richness < 3, b == 100) |>
  select(year, sample, island, cluster, poolname, 
         magna, longispina, pulex, nr_pools_within, 
         desiccation_dynamic) |>
  pivot_wider(names_from = sample, 
              values_from = c(magna, longispina, pulex)) |>
  remove_missing() |>
  mutate(isolation = if_else(nr_pools_within>median(nr_pools_within), "more surrounded", "less surrounded"),
         desicc = if_else(desiccation_dynamic>median(desiccation_dynamic), "high desiccation", "low desiccation")) |>
  summarise(mean = mean(magna_summer), 
            var = var(magna_summer),
            n = n(),
            .by = c(magna_spring, pulex_spring, year)) |>
  mutate(var_p = mean*(1-mean))

ggplot(test) + 
  aes(x = var, y = var_p) + 
  geom_point() +
  geom_abline(slope = 1)

pi <- 0.8
n  <- 100
dummy_data <- rbinom(n, 1, pi)
mean(dummy_data)
var(dummy_data)
pi*(1-pi)/n



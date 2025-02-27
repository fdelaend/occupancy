
## stats ----
## load the data
counts_env <- readRDS("../data/dataDaphnids.rds")
# variables are: 
# year: year of sampling
# sample: spring or summer sample
# island: unique island ID
# pool: unique pool ID
# cluster: unique cluster ID. A cluster is a collection of pools lying within a given distance from each other
# d: distance around the pool within which we are counting the number of surrounding pools
# richness: Daphnia species richness
# pools_dens: number of pools per running m of a cluster's perimeter
# desiccation_dynamic: estimates of pool desiccation rate
# latitude_corr and longitude_corr: lat and long
# nr_pools_within: nr of pools within a certain distance d

# Do some plotting: 
# seems like greater richness at higher nr of surr pools?
ggplot(counts_env) +
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = nr_pools_within, y = desiccation_dynamic, 
      col = as.factor(richness)) + 
  geom_point() + 
  facet_grid(d~sample, scales = "free", labeller = label_value) +
  labs(x="nr of surrounding rockpools", y = "desiccation rate", 
       col="richness")

# Regular GLM
model <- glm(richness ~ nr_pools_within + desiccation_dynamic, 
             data = counts_env |> filter(sample == "summer", d == 0.001), 
             family = poisson(link = "log"))
summary(model)

# Inspect residuals: probably spatial autocorrelation? 
counts_env |>
  mutate(pred = predict(model, newdata = counts_env, type = "response"),
         res = (pred - richness)^2) |>
  remove_missing() |>
  ggplot() +
  aes(x = latitude_corr, y = longitude_corr, colour = log10(res)) +
  geom_point() +
  facet_grid(.~sample, scales = "free", labeller = label_value)

# LEFTOVERS ----------

# Result 1: effect of mean desiccation: more single sp patches, fewer pairs, a bit more triplets
# Result 2: effect of nr of patches: fewer single sp patches 
# Interpretation: 
# Homogeneity of environmental conditions? Or pH, conductivity, plants not most important? Heterogeneity similar across all clusters
# Regional dominant? Not supported by the model? Pulex (and sometimes longi) according to the data?


# Hard to see so check out the slopes of fraction vs. predictor (p or desiccation)
slopes <- test |>
  filter(richness < 3) |>
  nest_by(sample, richness, year) %>%
  expand_grid(predictor = c("p", "desiccation_static_mean_mean")) %>%
  #mutate(n = map_dbl(data, ~ nrow(.x |> remove_missing()))) |>
  #filter(n > 0) |>
  mutate(formula = map(predictor, ~paste("fraction~",.x))) %>%
  mutate(slope = map2_dbl(data, formula, ~lm(data=.x, formula=as.formula(.y))$coefficients[[2]]))

ggplot(slopes) + 
  theme_bw() +
  aes(x=as.factor(richness), y=slope) + 
  geom_boxplot() + 
  facet_grid(predictor~sample, scales = "free", labeller = label_both) +
  geom_hline(yintercept = 0) +
  labs(x = "richness", y = "effect")
# ggsave("../figures/slope.pdf", width=6, height=4)

#visualize conductivity across systems; can do the same for pH
pH_conductivity |> 
  left_join(pools, by = join_by(island, pool)) |>
  filter(!is.na(latitude_corr)) |> #issue with island name resulting in coordinate not found
  ggplot() + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  scale_shape_manual(values = c(1:7)) +
  aes(x=latitude_corr, y=longitude_corr, col=log_conduct, 
      pch=cluster) + 
  geom_point() 

#Now for plant cover
plant_cover |>
  left_join(pools, by = join_by(island, pool)) |>
  filter(!is.na(latitude_corr)) |> #issue with island name resulting in coordinate not found
  ggplot() + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  scale_shape_manual(values = c(1:7)) +
  aes(x=latitude_corr, y=longitude_corr, col=plants, 
      pch=cluster) + 
  geom_point() 

# compute heterogenity
heterogeneity <- pH_conductivity |>
  #add cluster ID
  left_join(island_measures, by = join_by(island)) |>
  # environmental heterogeneity across pools within a cluster
  summarize(heterogeneity = sqrt(1/n()*sum((pH-mean(pH, na.rm=T))^2+
                                             (log_conduct-mean(log_conduct, na.rm=T))^2, na.rm=T)), 
            .by = cluster) |>
  remove_missing() |>
  mutate(heterogeneity_cut = factor(if_else(heterogeneity>1, "high", "low")))

plant_cover <- read_csv("../data/Ebert/Frederik_data1/plantcover_2013_2017.csv") |>
  remove_missing() |>
  # normalise across whole data set
  summarize(plants = mean(plants_rank), 
            .by = c("island", "pool")) |>
  mutate(plants = (plants-mean(plants))/sd(plants))
#plant data only for two years, also made temporal mean 

pH_conductivity <- read_csv("../data/Ebert/Frederik_data1/MetapopData_pH_conductivity_1998_2017.csv") |>
  # compute mean data across sampling occasions for every pool 
  summarize(pH = mean(pH, na.rm = T), 
            log_conduct = mean(log10(conduct_uS), na.rm = T), 
            .by = c("island", "pool")) |>
  # normalise across whole data set
  mutate(pH = (pH-mean(pH, na.rm = T))/sd(pH, na.rm = T),
         log_conduct = (log_conduct-mean(log_conduct, na.rm = T))/sd(log_conduct, na.rm = T)) 

#location of pools
pools           <- read_csv("../data/Ebert/Frederik_data1/Pools_coordinates_2017_vers7.csv") |>
  left_join(island_measures, by="island") |>
  select(island, pool, latitude_corr, longitude_corr, 
         cluster, number_pools_approx, island_perimeter_meter, pools_dens) 

# Fit the models 
# Resample 10 pools per cluster
resampled_data <- counts_env |> 
  group_by(cluster, year) |>
  nest() |>
  expand_grid(iteration = c(1:100)) |>
  mutate(data_sample = map(data, ~ .x[sample(1:nrow(.x), 12),])) |>
  select(!data) |>
  unnest(data_sample) |>
  group_by(iteration) |>
  nest() |>
  mutate(model = map(data, ~ tidy(glm(richness ~ pools_dens + desiccation_dynamic, 
                                      data = .x, 
                                      family = poisson(link = "log"))))) |>
  select(!data) |>
  unnest(model) |>
  filter(term != "(Intercept)")

ggplot(resampled_data) +
  aes(col = term, y = estimate) +
  geom_density()

#environmental variables
pH_conductivity <- read_csv("../data/Ebert/Frederik_data1/MetapopData_pH_conductivity_1998_2017.csv") |>
  # compute mean data across sampling occasions for every pool 
  summarize(pH = log10(mean(10^pH, na.rm = T)), 
            log_conduct = mean(log10(conduct_uS), na.rm = T), 
            .by = c("year", "island", "pool"))

# Check correlations among potential predictors of richness
counts_env |> 
  select(year, sample, nr_pools_within, desiccation_dynamic) |>
  cor(use = "pairwise.complete.obs")
# desiccation_dynamic and year seem to be correlated and 
# therefore I do not put them in the same model. 
# It is well known that pools dry more often in recent years,
# and so I assume that a year effect would mostly be a dryness effect.
# I further drop pH from the list of predictors because we don't have 
# pH measurements for 90% of the data. 
# Finally, I focus on the summer sample because most representative. 


# PLOT DATA 
ggplot(counts_env) +
  theme_bw() +
  #scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x = (nr_pools_within), 
      y = pairs) + 
  geom_boxplot() + 
  labs(x="nr of surrounding rockpools", 
       y = "species pairs") +
  facet_grid(b~sample, labeller = label_both)


pools |> 
  select(cluster, distances) |>
  unnest(distances) |>
  ggplot() + 
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = log10(distance), col = cluster) + 
  geom_density() +
  labs(x="log10(pairwise distance)", 
       y = "estimated density")

ggsave(filename = "../figures/distances.pdf", 
       width=4, height = 3, device = "pdf")
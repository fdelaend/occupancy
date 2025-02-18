
library(tidyverse)
library(broom)
library(lme4)

# READ DATA ----
cluster_measures <- read_csv("../data/Ebert/Frederik_data1/Island_measures_vers4.csv") |>
  mutate(cluster = factor(group_100)) |>
  mutate(nr_pools = sum(number_pools_approx),
         area = sum(island_area_sq_meter),
         peri = sum(island_perimeter_meter),
         .by = cluster) |>
  mutate(pools_dens = nr_pools / peri) |>
  select(cluster, island, pools_dens) 

pools <- read_csv("../data/Ebert/Frederik_data1/Pools_coordinates_2017_vers7.csv") |>
  filter(!grepl("A", pool), !grepl("A", name)) |> #island present twice in the data
  select(island, pool, latitude_corr, longitude_corr) |>
  left_join(cluster_measures |> select(cluster, island), 
            by = "island") |>
  group_by(cluster) |>
  nest() |>
  mutate(distances = map(data, ~ tidy(dist(.x |> select(latitude_corr, longitude_corr), upper = T)))) |>
  expand_grid(d = c(0.0001, 0.0005, 0.001, 0.002)) |>
  # nr of pools within 0.001 distance
  mutate(nr_pools_within = map2(distances, d, ~ .x |> 
                                 summarise(nr = sum((distance < .y)), .by = c(item1)))) |>
  # glue these results to the data
  mutate(data = map2(data, nr_pools_within, ~ .x |> mutate(nr_pools_within = .y$nr))) |>
  select(cluster, d, data) |>
  unnest(data)

counts          <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv") |>
  filter(!grepl("A", pool), !grepl("A", poolname)) |> #islands present twice in the data
  remove_missing() |>
  left_join(cluster_measures, by = join_by(island)) |>
  left_join(pools, by = join_by(cluster, island, pool))

desiccation        <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
  expand_grid(year = c(min(counts$year):max(counts$year))) |>
  left_join(read_csv("../data/Ebert/hydro.csv"), 
            by = join_by(poolname == pool, year)) |> 
  separate_wider_delim(poolname, delim = "-", names = c("island", "pool")) |>
  left_join(cluster_measures, by="island") |>
  summarise(desiccation_dynamic  = mean(predicted_desiccation_events, na.rm = T),
            desiccation_static_mean_mean   = mean(meandesi, na.rm = T),
            desiccation_static_mean_median = mean(mediandesi, na.rm = T),
            desiccation_static_median_median = median(mediandesi, na.rm = T),
            desiccation_static_median_mean = median(meandesi, na.rm = T),
            .by = c(year, cluster))

# add environmental data to count data
counts_env <- counts |>
  mutate(richness = magna + longispina + pulex) |>
  #lag one year to join with correct desiccation data in next step
  mutate(year_lag = year - 1) |>
  #add desiccation data
  left_join(desiccation, by = join_by(cluster, year_lag == year)) |>
  #remove na-s
  filter(!is.na(richness)) |>
  mutate(sample = if_else(sample == 1, "spring", "summer")) |> 
  filter(richness >0)

## save the data for olivia
counts_env |> 
  select(year, sample, island, pool, cluster, 
         richness, pools_dens, desiccation_dynamic, 
         latitude_corr, longitude_corr, nr_pools_within) |>
  saveRDS("../data/dataDaphnids.rds")

ggplot(counts_env) +
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = nr_pools_within, y = desiccation_dynamic, 
      col = as.factor(richness)) + 
  geom_point() + 
  facet_grid(d~sample, scales = "free", labeller = label_value) +
  labs(x="nr of surrounding rockpools", y = "desiccation rate", 
       col="richness")

## stats ----
## load the data
counts_env <- readRDS("../data/dataDaphnids.rds")
# variables are: 
# year: year of sampling
# sample: spring or summer sample
# island: unique island ID
# pool: unique pool ID
# cluster: unique cluster ID. A cluster is a collection of pools lying within a given distance from each other
# richness: Daphnia species richness
# pools_dens: number of pools per running m of a cluster's perimeter
# desiccation_dynamic: estimates of pool desiccation rate
# latitude_corr and longitude_corr: lat and long
# nr_pools_within: nr of pools within a certain distance

# Regular GLM
model <- glm(richness ~ nr_pools_within*sample + desiccation_dynamic*sample, 
             data = counts_env |> filter(d == 0.001, richness >0), 
             family = poisson(link = "log"))
summary(model)

# Inspect residuals
counts_env |>
  mutate(pred = predict(model, newdata = counts_env, type = "response"),
         res = (pred - richness)^2) |>
  ggplot() +
  aes(x = latitude_corr, y = longitude_corr, col = res) +
  geom_point() +
  facet_grid(.~sample, scales = "free", labeller = label_both)

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

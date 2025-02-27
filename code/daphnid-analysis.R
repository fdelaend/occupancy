
library(tidyverse)
library(broom)
library(lme4) # Mixed models
library(gstat) # variogram
library(sp) # coordinates
library(DHARMa) # glmer residuals check

# READ AND TREAT DATA ----
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
  expand_grid(d = c(1e-5, 1e-4, 0.0005, 0.001, 0.002)) |>
  # nr of pools within distance d
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
  remove_missing() |>
  mutate(sample = if_else(sample == 1, "spring", "summer"),
         monos = as.factor(if_else(richness==1, "present", "absent")),
         pairs = as.factor(if_else(richness==2, "present", "absent")),
         triplets = as.factor(if_else(richness==3, "present", "absent"))) |>
  filter(richness > 0, richness < 3) 
# rationale: study probability to observe 
# pairs vs monocultures,
# in pools that are viable. 
# triplets are extremely rare: 84 out of 135948 records

# PLOT DATA ------
ggplot(counts_env) +
  theme_bw() +
  #scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x = nr_pools_within) + 
  geom_histogram() + 
  labs(x="nr of surrounding rockpools") +
  facet_grid(d~sample, labeller = label_both)

# Species pairs seem more prevalent in pools that are well surrounded
# No effect of desiccation visible
ggplot(counts_env) +
  theme_bw() +
  #scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x = (nr_pools_within), 
      y = pairs) + 
  geom_boxplot() + 
  labs(x="nr of surrounding rockpools", 
       y = "species pairs") +
  facet_grid(d~sample, labeller = label_both)

## stats ----
# Fit models ----
counts_env_models <- counts_env |>
  filter(d>1e-5) |> #useless to try and fit this one
  group_by(sample, d) |>
  nest() |>
  mutate(model1 = map(data, ~ glmer(pairs ~ nr_pools_within + desiccation_dynamic +
                                      (1|island) + (1|year), data = .x, 
                                    family = binomial(link = "logit"))),
         model2 = map(data, ~ glmer(pairs ~ nr_pools_within + desiccation_dynamic +
                                      (1|island/pool) + (1|year), data = .x, 
                                    family = binomial(link = "logit"))),
         model_selection = map2(model1, model2, ~ tidy(anova(.x, .y, test = "Chisq"))),#, # model2 always better
         coefficients = map(model2, ~ coefficients(.x)))
         #dispersion = map(model2, ~ testDispersion(.x, plot = F)),
         #residuals = map(model2, ~ simulateResiduals(.x, use.u = T, plot = F))) #|> 

# Test residuals for some combo of d and sample
selected_model <- counts_env_models |> 
  filter(sample == "summer", 
         d == 1e-3) |> #1e-4, 5e-4, 1e-3, 2e-3
  ungroup() 

model <- selected_model$model2[[1]]
data <- selected_model$data[[1]]
testDispersion(model) 
sim <- simulateResiduals(fittedModel = model, 
                         plot = F, use.u = T)
plot(sim) 
plotResiduals(sim, data$year)
plotResiduals(sim, data$desiccation_dynamic)
plotResiduals(sim, data$nr_pools_within)
plotResiduals(sim, data$island)

data_sp <- data
coordinates(data) <- ~ latitude_corr + longitude_corr  # Define spatial coordinates
data$res <- residuals(sim)
variog <- variogram(res ~ 1, data)
plot(variog) #ok

summary(model)

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

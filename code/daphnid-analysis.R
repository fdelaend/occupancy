
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
  mutate(data = map2(data, nr_pools_within, ~ .x |> mutate(nr_pools_within = .y$nr))) 

pools |> 
  select(cluster, distances) |>
  unnest(distances) |>
  ggplot() + 
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = log10(distance), col = cluster) + 
  geom_density() +
  labs(x="pairwise distance", 
       y = "estimated density")
  
ggsave(filename = "../figures/distances.pdf", 
       width=4, height = 3, device = "pdf")

pools <- pools |>
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

## stats ----
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
         model3 = map(data, ~ glm(pairs ~ nr_pools_within, data = .x, 
                                    family = binomial(link = "logit"))),
         model_selection = map2(model1, model2, ~ tidy(anova(.x, .y, test = "Chisq"))),#, # model2 always better
         coefficients = map(model2, ~ coefficients(.x)))
         #dispersion = map(model2, ~ testDispersion(.x, plot = F)),
         #residuals = map(model2, ~ simulateResiduals(.x, use.u = T, plot = F))) #|> 

# Test residuals for some combo of d and sample
selected_model <- counts_env_models |> 
  filter(sample == "summer", 
         d == 5e-4) |> #1e-4, 5e-4, 1e-3, 2e-3
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

##plot ----
counts_env_plot <- counts_env_models |>
  mutate(data = map2(data, model3, ~ .x |> 
                       mutate(preds = predict(.y, type="response")))) |>
  select(!contains("model") & !coefficients) |>
  unnest(data) |>
  ungroup() |>
  summarise(prop = mean(pairs=="present"),
            pred = mean(preds),
            pred_sd = sd(preds),
            p = mean(nr_pools_within),
            desiccation_dynamic = mean(desiccation_dynamic),
            .by = c(year, sample, cluster, d))

ggplot(counts_env_plot) +
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = p, y = prop, 
      col = as.factor(year)) +
  geom_point(show.legend = FALSE) + 
  geom_line(aes(x = p, y = pred), col = "black", show.legend = FALSE) +
  labs(x="nr of surrounding rockpools", 
       y = "proportion of pools with pairs") +
  facet_grid(d~sample, labeller = label_both) 

ggsave(filename = "../figures/glm.pdf", 
       width=5, height = 5, device = "pdf")

#leftovers
# PLOT DATA ------
ggplot(counts_env) +
  theme_bw() +
  #scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x = (nr_pools_within), 
      y = pairs) + 
  geom_boxplot() + 
  labs(x="nr of surrounding rockpools", 
       y = "species pairs") +
  facet_grid(d~sample, labeller = label_both)
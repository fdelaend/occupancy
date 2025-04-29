library(grateful)
library(tidyverse)
library(xtable)
library(broom.mixed)
library(lme4) # Mixed models
library(gstat) # variogram
library(geodist) #distance calculation
library(sp) # coordinates
library(DHARMa) # glmer residuals check

# READ AND TREAT DATA ----
cluster_measures <- read_csv("../data/Ebert/Frederik_data1/Island_measures_vers4.csv") |>
  mutate(cluster = factor(group_100)) |>
  select(cluster, island) 

pools <- read_csv("../data/Ebert/Frederik_data1/Pools_coordinates_2017_vers7.csv") |>
  filter(!grepl("A", pool), !grepl("A", name)) |> #island present twice in the data
  select(island, pool, latitude_corr, longitude_corr) |>
  left_join(cluster_measures |> select(cluster, island), 
            by = "island") |>
  group_by(cluster) |>
  nest() |>
  mutate(distances = map(data, ~ geodist(.x |> select(latitude_corr, longitude_corr)))) |>
  expand_grid(b = c(10, 50, 100, 200)) |>
  # nr of pools within distance b
  mutate(nr_pools_within = map2(distances, b, ~ rowSums(.x < .y)-1)) |>
  # glue these results to the data
  mutate(data = map2(data, nr_pools_within, ~ .x |> mutate(nr_pools_within = .y))) 

pools <- pools |>
  select(cluster, b, data) |>
  unnest(data)

counts          <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv") |>
  filter(!grepl("A", pool), !grepl("A", poolname)) |> #islands present twice in the data
  remove_missing() |>
  left_join(cluster_measures, by = join_by(island)) |>
  left_join(pools, by = join_by(cluster, island, pool), 
            relationship = "many-to-many")

desiccation        <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
  expand_grid(year = c(min(counts$year):max(counts$year))) |>
  left_join(read_csv("../data/Ebert/hydro.csv"), 
            by = join_by(poolname == pool, year)) |> 
  separate_wider_delim(poolname, delim = "-", names = c("island", "pool")) |>
  left_join(cluster_measures, by="island") |>
  summarise(desiccation_dynamic  = mean(predicted_desiccation_events, na.rm = T),
            .by = c(year, cluster))

# add desiccation data to count data
counts_des <- counts |>
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
         triplets = as.factor(if_else(richness==3, "present", "absent"))) 

## stats ----
# rationale: study probability to observe 
# pairs vs monocultures,
# in pools that are viable. 
# Remove triplets bc extremely rare: 21 out of 6297 occurrences of non-empty pools (0.33%):
counts_des |>
  filter(b == 10, richness == 3) |> #richness == 3, or >0
  nrow()

# Fit glmer for different values of b
counts_des_models <- counts_des |>
  filter(richness > 0, richness < 3) |>
  group_by(sample, b) |>
  nest() |>
  mutate(model1 = map(data, ~ glmer(pairs ~ nr_pools_within + desiccation_dynamic +
                                      (1|island) + (1|year), data = .x, 
                                    family = binomial(link = "logit"))),
         model2 = map(data, ~ glmer(pairs ~ nr_pools_within + desiccation_dynamic +
                                      (1|island/pool) + (1|year), data = .x, 
                                    family = binomial(link = "logit"))),
         # Simplified model used for plotting
         model3 = map(data, ~ glm(pairs ~ nr_pools_within, data = .x, 
                                    family = binomial(link = "logit"))),
         model_selection = map2(model1, model2, ~ tidy(anova(.x, .y, test = "Chisq"))))#, # model2 always better
         #coefficients = map(model2, ~ coefficients(.x)))
         #dispersion = map(model2, ~ testDispersion(.x, plot = F)),
         #residuals = map(model2, ~ simulateResiduals(.x, use.u = T, plot = F))) #|> 

###Model selection table and dispersal test model 2----
counts_des_models |>
  ungroup() |>
  #do dispersal test of model 2
  mutate(`p, disp.` = map(model2, ~ testDispersion(.x, plot = F)$p.value)) |>
  select(!data & !model1 & !model2 & !model3) |>
  unnest(model_selection) |>
  mutate(model = if_else(term == ".x", "model 1", "model 2")) |>
  select(!term & !npar) |>
  mutate(b = as.character(b),
         `p, comp.` = as.character(signif(p.value,3))) |>
  select(-p.value) |>
  relocate(model, .after = sample) |>
  relocate(`p, disp.`, .after = deviance) |>
  as.data.frame() |>
  xtable(digits = c(0, 0, 0, 2, rep(0, 4), 3, 0, 2, 2)) |>
  print(type = "latex",
        scientific = TRUE,
        include.rownames = FALSE)

###Residual analysis ----

for (i in 1:nrow(counts_des_models)) {
  model <- counts_des_models$model2[[i]]
  data <- counts_des_models$data[[i]] |>
    rename(desiccation = desiccation_dynamic,
           `nr of pools` = nr_pools_within)
  
  sim <- simulateResiduals(fittedModel = model, 
                           plot = F, use.u = T)
  
  pdf(paste0("../figures/", i, "-resid-qq.pdf"), 
      width = 8, height = 5)
  plot(sim) 
  dev.off()
  
  pdf(paste0("../figures/", i, "-resid.pdf"), 
      width = 10, height = 9)
  par(mfrow = c(2, 3))
  plotResiduals(sim, as.factor(data$year))
  title(xlab = "catPred", col.lab = "white")
  title(xlab = "year")
  plotResiduals(sim, data$desiccation)
  plotResiduals(sim, data$`nr of pools`)
  plotResiduals(sim, as.factor(data$island))
  title(xlab = "catPred", col.lab = "white")
  title(xlab = "island")
  coordinates(data) <- ~ latitude_corr + longitude_corr  # Define spatial coordinates
  data$res <- residuals(sim)
  variog <- variogram(res ~ 1, data)
  plot(variog$dist, variog$gamma, main = "variogram",
       xlab = "distance", ylab = "semivariance", 
       ylim = c(0, max(variog$gamma)))
  dev.off()
}

###estimated effects of final models ----
counts_des_models |>
  ungroup() |>
  select(!data & !model1 & !model3 & !model_selection) |>
  mutate(results = map(model2, ~ tidy(.x))) |>
  select(!model2) |>
  unnest(results) |>
  select(!effect) |>
  mutate(b = as.character(b),
         estimate = as.character(signif(estimate,3)),
         std.error = as.character(signif(std.error,3)),
         statistic = as.character(signif(statistic,3)),
         p.value = as.character(signif(p.value,3))) |>
  as.data.frame() |>
  xtable() |>
  print(type = "latex",
        scientific = TRUE,
        include.rownames = FALSE)

##plot ----
counts_des_plot <- counts_des_models |>
  mutate(data = map2(data, model3, ~ .x |> 
                       mutate(preds = predict(.y, type="response")))) |>
  select(!contains("model")) |>
  unnest(data) |>
  ungroup() |>
  summarise(prop = mean(pairs=="present"),
            pred = mean(preds),
            pred_sd = sd(preds),
            p = mean(nr_pools_within),
            n = length(nr_pools_within), # Weight for the logistic regression on proportions
            desiccation_dynamic = mean(desiccation_dynamic),
            .by = c(year, sample, cluster, b))

ggplot(counts_des_plot |> filter(b == 100)) +
  theme_bw() +
  scale_x_reverse() + #for talk
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = p, y = prop, 
      col = as.factor(b),
      weight = n) + # Weight for the logistic regression on proportions
  geom_point(alpha = 0.2, show.legend = FALSE) + 
  geom_point(shape = 1, alpha = 0.5, show.legend = FALSE) + 
  geom_smooth(aes(group = as.factor(b)), # White border to highlight lines
              method = "glm", 
              method.args = list(family = binomial(link = "logit")),
              se = F, 
              color = "white", 
              lwd = 1.4,
              lineend='round') +  # Round edge of lines
  geom_smooth(method = "glm", 
              method.args = list(family = binomial(link = "logit")),
              se = F,
              lineend='round',
              lwd = 0.8) +
  labs(x="mean nr of pools within distance, by cluster", 
       y = "proportion of cluster pools with pairs", 
       col = "distance (m)") +
  facet_grid(.~sample, labeller = label_both) 

ggsave(filename = "../figures/glm-talk.pdf", 
       width=5, height = 3, device = "pdf")



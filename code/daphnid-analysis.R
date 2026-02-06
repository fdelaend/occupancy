library(grateful)
library(xtable)
library(broom.mixed)
library(lme4) # Mixed models
library(gstat) # variogram
library(geodist) #distance calculation
library(sp) # coordinates
library(DHARMa) # glmer residuals check
library(ggeffects) # plot effects of GLMMs
source("tools-other.R") # function to fit glmms for priority effects

# Read and organise data ------------

## Cluster data: List of islands and clusters they belong to -------
cluster_measures <- read_csv("../daphnid-data/island-measures.csv") |>
  mutate(cluster = factor(group_100)) |>
  select(cluster, island) 

## Spatial location of pools and nr of surrounding pools within cluster --------- 
pool_coords <- read_csv("../daphnid-data/pool-coords.csv") |>
  # Pool coordinates 
  filter(!grepl("A", pool), !grepl("A", name)) |> #island present twice in the data
  select(island, pool, name, latitude_corr, longitude_corr) |>
  # Add cluster info 
  left_join(cluster_measures |> select(cluster, island), 
            by = "island") |>
  group_by(cluster) |>
  nest() |>
  # Compute nr of pools within a distance "b"
  mutate(distances = map(data, ~ geodist(.x |> select(latitude_corr, longitude_corr)))) |>
  expand_grid(b = c(50, 100)) |>
  # nr of pools within distance b
  mutate(nr_pools_within = map2(distances, b, ~ rowSums(.x < .y)-1)) |>
  # glue these results to the data
  mutate(data = map2(data, nr_pools_within, ~ .x |> mutate(nr_pools_within = .y))) |>
  # keep only what we need
  select(cluster, b, data) |>
  # unnest
  unnest(data)

## Daphnid data -------
# Read Daphnid data and join cluster info and nr of pools within "b"
daphnids          <- read_csv("../daphnid-data/daphnid-counts.csv") |>
  filter(!grepl("A", pool), !grepl("A", poolname)) |> #islands present twice in the data
  drop_na() |>
  left_join(cluster_measures, by = join_by(island)) |>
  left_join(pool_coords, by = join_by(cluster, island, pool), 
            relationship = "many-to-many") |>
  #compute daphnid richness
  mutate(richness = magna + longispina + pulex)

## Desiccation data ---------
desiccation        <- tibble(name = unique(pool_coords$name)) |>
  expand_grid(year = c(min(daphnids$year):max(daphnids$year))) |>
  left_join(read_csv("../daphnid-data/desiccation.csv"), 
            by = join_by(name == pool, year)) |> 
  separate_wider_delim(name, delim = "-", names = c("island", "pool")) |>
  left_join(cluster_measures, by="island") |>
  # desiccation is a cluster-level predictor, computed as the mean across pools within a cluster
  summarise(desiccation_dynamic  = mean(predicted_desiccation_events, na.rm = T),
            .by = c(year, cluster))

## Add desiccation data to daphnids data --------
daphnids_des <- daphnids |>
  #lag one year to join with correct desiccation data in next step
  mutate(year_lag = year - 1) |>
  #add desiccation data
  left_join(desiccation, by = join_by(cluster, year_lag == year)) |>
  #remove na-s
  drop_na() |>
  mutate(sample = if_else(sample == 1, "spring", "summer"),
         monos = as.factor(if_else(richness==1, "present", "absent")),
         pairs = as.factor(if_else(richness==2, "present", "absent")),
         triplets = as.factor(if_else(richness==3, "present", "absent"))) 

# Statistics ----------
# Remove triplets bc extremely rare: 21 out of 6297 occurrences of non-empty pools (0.33%):
daphnids_des |>
  filter(b == 100, richness == 3) |> #richness == 3, or >0
  nrow()

daphnids_des <- daphnids_des |>
  filter(richness > 0, richness < 3)

## Priority effects ------
### Prep data ------
data_priority <- daphnids_des |>
  # filter data: only isolated pools in which a single species 
  # was present in spring. Pick some value of b, not important, because not used. 
  filter(!(sample == "spring" & richness != 1),
         b == 100, nr_pools_within < 20) |>
  select(year, sample, island, cluster, pool, latitude_corr, longitude_corr,
         magna, longispina, pulex) |>
  # create variables that code for presence of a sp in a season
  pivot_wider(names_from = sample, 
              values_from = c(magna, longispina, pulex)) |>
  drop_na() |>
  #create variables that code for presence of a sp in summer w/o other sp being present
  mutate(magna_summer_alone      = magna_summer * (1-longispina_summer) * (1-pulex_summer),
         longispina_summer_alone = (1-magna_summer) * longispina_summer * (1-pulex_summer),
         pulex_summer_alone      = (1-magna_summer) * (1-longispina_summer) * pulex_summer) |>
  #convert to factors for stats
  mutate(island = as.factor(island),
         year = as.factor(year))

all_species <- c("magna", "longispina", "pulex")

### Model fitting --------
models_priority <- tibble(focal = all_species) |>
  mutate(model = map(focal, 
                     ~ fit_priority_2(focal = .x, adapt_delta = 0.999,
                                      data = data_priority))) 

saveRDS(models_priority, "../data/models_priority.RDS")
models_priority <- readRDS("../data/models_priority.RDS")

### Model plotting --------
# get out the typical names of the parameters we want to plot
colnames(as.matrix(models_priority$model[[1]]))
# do the plotting
models_priority |>
  mutate(draws = map(model, ~as_draws_df(.x) |>
                       select(1:5) |>
                       rename(intercept = 1,
                              `alone in spring` = 2,
                              `sd island` = 3,
                              `sd island:pool` = 4,
                              `sd year` = 5))) |> #ignore warnings
  select(!model) |>
  unnest(draws) |>
  pivot_longer(2:6, names_to = "parameter", values_to = "value") |>
  ggplot() + 
  aes(x = value, y = focal, fill = parameter) + 
  stat_halfeye(alpha = 0.6, slab_color = NA) +
  #facet_wrap(~ parameter, scales = "free_x") +
  theme_bw()
  
ggsave(filename = "../figures/priority.pdf", 
       width=5, height = 3, device = "pdf")
#comparison w out of the box solution
#mcmc_areas(models_priority$model[[3]], 
#           pars = colnames(as.matrix(models_priority$model[[3]]))[1:5])

### Model checking --------
test <- models_priority |>
  mutate(ppc_draws = map(model, ~posterior_predict(.x)))

p1 <- ppc_dens_overlay(y = data_priority$magna_summer_alone, 
                 yrep = test$ppc_draws[[1]][1:100, ]) + 
  labs(title = "magna") 
p2 <- ppc_dens_overlay(y = data_priority$longispina_summer_alone, 
                 yrep = test$ppc_draws[[2]][1:100, ]) +
  labs(title = "longispina") 
p3 <- ppc_dens_overlay(y = data_priority$pulex_summer_alone, 
                 yrep = test$ppc_draws[[3]][1:100, ]) +
  labs(title = "pulex")

(p1 | p2 | p3) 

ggsave("../figures/check_priority.pdf",
       width=8, height = 3, device = "pdf")                 
                               
## Proportion of pairs -----
### Model fitting ------
models_prop <- daphnids_des |>
  group_by(sample, b) |>
  nest() |>
  mutate(model = map(data, ~ brm(pairs ~ nr_pools_within + desiccation_dynamic +
                                      (1|island/pool) + (1|year), data = .x, 
                                 family = brms::bernoulli(),
                                 control = list(adapt_delta = 0.999),
                                 warmup = 6000,
                                 iter = 8000)))

saveRDS(models_prop, "../data/models_prop.RDS")
models_prop <- readRDS("../data/models_prop.RDS")

### Model plotting --------
#### Posteriors -------
# get out the typical names of the parameters we want to plot
colnames(as.matrix(models_prop$model[[1]]))
# do the plotting
models_prop |>
  mutate(draws = map(model, ~as_draws_df(.x) |>
                       select(1:6) |>
                       rename(intercept = 1,
                              `nr of nearby pools` = 2,
                              `desiccation` = 3,
                              `sd island` = 4,
                              `sd island:pool` = 5,
                              `sd year` = 6))) |> #ignore warnings
  select(!model & !data) |>
  unnest(draws) |>
  pivot_longer(3:8, names_to = "parameter", values_to = "value") |>
  ggplot() + 
  aes(x = value, y = as.factor(b), fill = parameter) + 
  stat_halfeye(alpha = 0.6, slab_color = NA) +
  facet_grid(sample~parameter, 
             scales = "free") +
  theme_bw() +
  labs(y = "b")

ggsave(filename = "../figures/proportion.pdf", 
       width=12, height = 3, device = "pdf")

#### Fixed effect -------
# prediction grid across nr_pools_within
new_data <- tibble(nr_pools_within = seq(1, 100,
    length.out = 200), 
    desiccation_dynamic = mean(daphnids_des$desiccation_dynamic, na.rm = TRUE))

# summary of observed data
daphnids_des_summary <- daphnids_des |>
  summarise(p_mean = mean(richness==2, na.rm=T), 
            nr_pools_within = mean(nr_pools_within),
            .by = c(year, sample, b, cluster))

test <- models_prop |>
  select(!data) |>
  # population-level predicted probabilities
  mutate(ep = map(model, ~ posterior_epred(.x, newdata = new_data, 
                                           re_formula = NA))) |>
  mutate(
    pred = map(ep, \(m) {
      # m is draws x grid
      tibble(
        nr_pools_within = new_data$nr_pools_within,
        desiccation_dynamic = new_data$desiccation_dynamic,
        p_mean = colMeans(m),
        p_lo   = apply(m, 2, quantile, probs = 0.025),
        p_hi   = apply(m, 2, quantile, probs = 0.975)
      )
    })
  ) |>
  select(-model, -ep) |>
  unnest(pred)

ggplot(test, aes(x = nr_pools_within, y = p_mean, col = as.factor(b))) +
  scale_color_viridis_d(option="plasma", end=0.9) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi, col = as.factor(b)), 
              alpha = 0.25, lwd=0.3) +
  geom_line() +
  facet_wrap(.~ sample) +
  labs(x = "nr of nearby pools", y = "Predicted probability", 
       col = "distance (m)") +
  theme_bw()

ggsave(filename = "../figures/proportion_nearby.pdf", 
       width=4, height = 2, device = "pdf")

### Model checking --------
test <- models_prop |>
  mutate(ppc_draws = map(model, ~posterior_predict(.x)))

for (b in c(50, 100)){
  for (sample in c("spring", "summer")) {
    index <- which((test$b==b)&(test$sample==sample))
    p <- ppc_dens_overlay(y = (test$data[[index]]$pairs=="present")*1, 
                          yrep = test$ppc_draws[[index]][1:100, ]) + 
      labs(title = paste("b=", b, "; sample = ", sample))
    assign(paste0("b",b,sample), p)
  }
}


(b50spring | b50summer)/(b100spring | b100summer) 

ggsave("../figures/check_prop.pdf",
       width=6, height = 5, device = "pdf")        

# LEFTOVERS --------------------
### Dispersal test ----
daphnids_des_models |>
  ungroup() |>
  #do dispersal test
  mutate(`p, disp.` = map(model, ~ testDispersion(.x, plot = F)$p.value)) |>
  select(-model) |>
  as.data.frame() |>
  xtable(digits = c(rep(0, 2), 3)) |>
  print(type = "latex",
        scientific = TRUE,
        include.rownames = FALSE)

### Residual analysis ------

coordinates(data) <- ~ latitude_corr + longitude_corr  # Define spatial coordinates

for (i in 1:nrow(daphnids_des_models)) {
  model <- daphnids_des_models$model[[i]]
  sim <- simulateResiduals(fittedModel = model, 
                           plot = F, use.u = T)
  
  pdf(paste0("../figures/", i, "-resid-qq-prior.pdf"), 
      width = 8, height = 5)
  plot(sim) 
  dev.off()
  
  pdf(paste0("../figures/", i, "-resid-prior.pdf"), 
      width = 8, height = 5)
  par(mfrow = c(1, 2))
  plotResiduals(sim, as.factor(data$year))
  title(xlab = "catPred", col.lab = "white")
  title(xlab = "year")
  data$res <- residuals(sim)
  variog <- variogram(res ~ 1, data)
  plot(variog$dist, variog$gamma, main = "variogram",
       xlab = "distance", ylab = "semivariance", 
       ylim = c(0, max(variog$gamma)))
  dev.off()
}

### Estimated effects of final models ----
daphnids_des_models |>
  ungroup() |>
  mutate(results = map(model, ~ tidy(.x))) |>
  select(!model) |>
  unnest(results) |>
  #select(!effect) |>
  mutate(estimate = as.character(signif(estimate,3)),
         std.error = as.character(signif(std.error,3)),
         statistic = as.character(signif(statistic,3)),
         p.value = as.character(signif(p.value,3))) |>
  as.data.frame() |>
  xtable() |>
  print(type = "latex",
        scientific = TRUE,
        include.rownames = FALSE)

### Plot -----
daphnids_des_models_plot <- daphnids_des_models |>
  mutate(predictors = map(model, ~ gsub("as.factor\\((.*)\\)", "\\1", attr(terms(.x), "term.labels"))),
         predictions = map2(model, predictors, ~ ggpredict(.x, terms = .y, bias_correction = FALSE))) |>
  select(!predictors) |>
  unnest(predictions) #|>
  #mutate(group = factor(group, levels = c(0, 1),
#                     labels = c("no", "yes")))

ggplot(daphnids_des_models_plot, 
       aes(x = factor(x), y = predicted)) + #,
  #colour = group, group = group
  #scale_colour_manual(values = c("grey50", "black")) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  #geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.1,
                position = position_dodge(width = 0.3)) +
  scale_y_continuous("P(focal sp. alone in summer)", 
                     limits = c(0, 1)) +
  scale_x_discrete("Focal sp. alone in spring", 
                   labels = c("0" = "no", "1" = "yes")) +
  #scale_colour_discrete("Other sp. in spring") +
  theme_bw() + 
  facet_grid(.~focal) #+
  #labs(col = "Other sp. in spring")

ggsave(filename = "../figures/priority.pdf", 
       width=5, height = 2, device = "pdf")

## Alternative option: -----
# focus only on isolates ponds, 
# that have a single sp in spring and either 
# another sp (exclusion) or both sp (coex) in summer. 
# Model the occurrence of both events.

daphnids_des_alt <- daphnids_des |>
  mutate(keep = if_else(((richness==1)&(sample=="spring"))|sample=="summer", 1, 0)) |>
  filter(richness > 0, richness < 3, keep == 1,
         b == 100, nr_pools_within < 20) |>
  select(richness, year, sample, island, cluster, poolname, latitude_corr, longitude_corr,
         magna, longispina, pulex) 

data <- daphnids_des_alt |>
  pivot_wider(names_from = sample, 
              values_from = c(magna, longispina, pulex)) |>
  drop_na() |>
  mutate(obs = row_number()) |>
  #filter(n>20) |>
  mutate(island = as.factor(island),
         year = as.factor(year))

all_species <- c("magna", "longispina", "pulex")





### Model selection table and dispersal test model 2----
daphnids_des_models |>
  ungroup() |>
  #do dispersal test of model 2
  mutate(`p, disp.` = map(version2, ~ testDispersion(.x, plot = F)$p.value)) |>
  select(!data & !version1 & !version2 & !version3) |>
  unnest(model_selection) |>
  mutate(model = if_else(term == ".x", "version 1", "version 2")) |>
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

### Residual analysis ----

for (i in 1:nrow(daphnids_des_models)) {
  model <- daphnids_des_models$model2[[i]]
  data <- daphnids_des_models$data[[i]] |>
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

### Estimated effects of final models ----
daphnids_des_models |>
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

### Plot ----
daphnids_des_plot <- daphnids_des_models |>
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

ggplot(daphnids_des_plot) +
  theme_bw() +
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

ggsave(filename = "../figures/glm.pdf", 
       width=5, height = 3, device = "pdf")
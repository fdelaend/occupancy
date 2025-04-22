library(grateful)
library(tidyverse)
library(xtable)
library(broom.mixed)
library(lme4) # Mixed models
library(gstat) # variogram
library(sp) # coordinates
library(DHARMa) # glmer residuals check

readRDS("Data.rds")

## stats ----
# rationale: study probability to observe 
# pairs vs monocultures,
# in pools that are viable. 
# Remove triplets bc extremely rare: 21 out of 6297 occurrences of non-empty pools (0.33%):
counts_env |>
  filter(b == 1e-5, richness == 3) |> #richness == 3, or >0
  nrow()

# Fit glmer for different values of b
counts_env_models <- counts_env |>
  filter(richness > 0, richness < 3) |>
  filter(b > 1e-5) |> #useless to try and fit this one
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
counts_env_models |>
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

for (i in 1:nrow(counts_env_models)) {
  model <- counts_env_models$model2[[i]]
  data <- counts_env_models$data[[i]] |>
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
counts_env_models |>
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
counts_env_plot <- counts_env_models |>
  mutate(data = map2(data, model3, ~ .x |> 
                       mutate(preds = predict(.y, type="response")))) |>
  select(!contains("model")) |>
  unnest(data) |>
  ungroup() |>
  summarise(prop = mean(pairs=="present"),
            pred = mean(preds),
            pred_sd = sd(preds),
            p = mean(nr_pools_within),
            desiccation_dynamic = mean(desiccation_dynamic),
            .by = c(year, sample, cluster, b))

ggplot(counts_env_plot) +
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = p, y = prop, 
      col = as.factor(b)) +
  geom_point(alpha = 0.5, show.legend = FALSE) + 
  geom_line(aes(x = p, y = pred, 
                col = as.factor(b)),
            lwd = 1.2) +
  labs(x="mean nr of surrounding rockpools within distance", 
       y = "proportion of pools with pairs", 
       col = "distance") +
  facet_grid(.~sample, labeller = label_both) 

ggsave(filename = "../figures/glm.pdf", 
       width=5, height = 3, device = "pdf")

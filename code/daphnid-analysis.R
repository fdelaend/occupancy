
library(tidyverse)
library(lme4)

# READ DATA ----
island_measures <- read_csv("../data/Ebert/Frederik_data1/Island_measures_vers4.csv") |>
  select(island, group_100, group_200) |>
  mutate(cluster = factor(group_100)) |>
  select(island, cluster)

counts          <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv") |>
  left_join(island_measures, by="island")

#desiccation     <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
desiccation_static <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
  expand_grid(year = c(min(counts$year):max(counts$year)))

desiccation        <- desiccation_static |>
  left_join(read_csv("../data/Ebert/hydro.csv"), 
            by = join_by(poolname == pool, year)) |> 
  separate_wider_delim(poolname, delim = "-", names = c("island", "pool")) |>
  left_join(island_measures, by="island") |>
  summarise(desiccation_dynamic  = mean(predicted_desiccation_events, na.rm = T),
            desiccation_static_mean_mean   = mean(meandesi, na.rm = T),
            desiccation_static_mean_median = mean(mediandesi, na.rm = T),
            desiccation_static_median_median = median(mediandesi, na.rm = T),
            desiccation_static_median_mean = median(meandesi, na.rm = T),
            .by = c(year, cluster))

ggplot(desiccation) +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x = desiccation_static_mean_mean, y = desiccation_dynamic, 
      col = as.factor(year)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

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
  select(island, pool, latitude_corr, longitude_corr, cluster) 

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

# PLOT DATA ----
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

# add desiccation to data and 
# compute richness and fraction
fractions <- counts |>
  mutate(richness = magna + longispina + pulex) |>
  #remove empty pools
  filter(richness > 0) |>
  #total nr of pools per cluster, across all sampling events
  mutate(p = length(unique(poolname)), .by = cluster) |>
  #number of pools evaluated in a given cluster at a given sampling event in a given year
  mutate(n = length(unique(poolname)), 
         .by = c(cluster, sample, year)) |> 
  #remove all cluster - sampling event - year combos that have too small n
  filter(n > 20) |>
  #the weight of a single pool is 1/n. 
  #Making the sum across richness levels seems a tidy way to get the fraction of pools per richness level
  mutate(weight = 1/n) |>
  #fraction of pools with given richness level
  summarise(fraction = sum(weight), 
            .by = c(cluster, p, sample, year, richness, n)) |> #, "author"
  #add desiccation data
  left_join(desiccation, by = join_by(cluster, year)) |>
  #remove na-s
  filter(!is.na(desiccation_dynamic))

#plot of fraction vs. nr of patches p or desiccation
fractions |>
  filter(richness < 3) |>
  mutate(sample = if_else(sample == 1, "spring", "summer")) |>
  ggplot() +
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=desiccation_dynamic, y=fraction, col=year)) + 
  facet_grid(richness~sample, scales = "free", 
             labeller = label_bquote(rows=paste("richness = ", .(richness)),
                                     cols=paste(.(sample), " sample"))) +
  aes(x=desiccation_dynamic, y=fraction, col=year, group=as.factor(year)) +
  geom_smooth(lwd=0.5, method = lm, se=F, show.legend = F) +
  labs(x = "nr of desiccation events (proxy for dispersal rate)", 
       y = "fraction")

ggsave(paste0("../figures/f_desiccation.pdf"), width=4.5, height = 3, 
       device = "pdf")

# STATS ----
# prep data for stats
data <- test |>
  filter(richness == 2, sample == 2) |>
  select(p, year, n, fraction, desiccation_dynamic) |>
  mutate(year = as.factor(year))

model <- glmer(data = data, 
               formula = fraction ~ desiccation_dynamic + (1|year), 
               family = binomial, weights = n)

summary(model)
plot(model)

#HERE ---- (double check this; why fewer data than for fractions??)
#Dominance of a given species?
dom <- counts |>
  mutate(richness = magna + longispina + pulex) |>
  #remove empty pools
  filter(richness > 0) |>
  #total nr of pools per cluster, across all sampling events
  mutate(p = length(unique(poolname)), .by = cluster) |>
  #number of pools evaluated in a given cluster at a given sampling event in a given year
  mutate(n = length(unique(poolname)), 
         .by = c(cluster, sample, year)) |> 
  #remove all cluster - sampling event - year combos that have too small n
  filter(n > 20) |>
  mutate(`only magna` = (magna)*(1-longispina)*(1-pulex),
         `only longi` = (1-magna)*(longispina)*(1-pulex),
         `only pulex` = (1-magna)*(1-longispina)*(pulex)) |>
  select(!magna & !longispina & !pulex & !water & !author & !richness) |>
  pivot_longer(starts_with("only"), 
               names_to = "species", values_to = "value") |>
  #fraction of pools with only species x
  summarize(fraction = sum(value) / mean(n), 
            .by = c(cluster, p, sample, year, species, n)) |> #, "author"
  #add desiccation data
  left_join(desiccation, by = join_by(cluster, year)) |>
  filter(!is.na(desiccation_dynamic))


dom |>
  mutate(sample = if_else(sample == 1, "spring", "summer")) |>
  ggplot() +
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=desiccation_dynamic, y=fraction, col=year)) +
  facet_grid(species~sample, scales = "free", 
             labeller = label_bquote(rows=paste(.(species)),
                                     cols=paste(.(sample), " sample"))) +
  aes(x=desiccation_dynamic, y=fraction, col=year, group=as.factor(year)) +
  geom_smooth(lwd=0.5, method = lm, se=F, show.legend = F) +
  labs(x = "nr of desiccation events (proxy for dispersal rate)", 
       y = "fraction")

ggsave(paste0("../figures/f_dom.pdf"), width=4.5, height = 4, 
       device = "pdf")

ggplot(dom) + 
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x=desiccation, y=fraction, col=as.factor(year)) + 
  geom_point(show.legend = F) + 
  facet_grid(sample~species) +
  geom_smooth(lwd=0.5, method = lm, se=F, show.legend = F)
# ggsave("../figures/dom.pdf", width=6, height=4)

data <- dom |> filter(species == "only_longi", sample == 2) |> 
  select(fraction, n, desiccation_dynamic, year)
# try a model 
model <- glmer(fraction ~ desiccation_dynamic + 
                 (1 + desiccation_dynamic | year), 
               data = data, weight = n, 
               family = binomial(link = "logit"))
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


library(tidyverse)

# READ DATA ----
island_measures <- read_csv("../data/Ebert/Frederik_data1/Island_measures_vers4.csv") |>
  select(island, group_100, group_200) |>
  mutate(cluster = factor(group_100)) |>
  select(island, cluster)

counts          <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv") |>
  left_join(island_measures, by="island")

#desiccation     <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
desiccation     <- read_csv("../data/Ebert/hydro.csv") |>
  separate_wider_delim(pool, delim = "-", names = c("island", "pool")) |>
  left_join(island_measures, by="island") |>
  summarise(desiccation = mean(predicted_desiccation_events, na.rm = T), 
            .by = c(year, cluster)) |>
  remove_missing()

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

# ANALYSE DATA ----
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

# add heterogeneity and desiccation to data and 
# compute richness and fraction
test<-counts |>
  mutate(richness = magna + longispina + pulex) |>
  #remove empty pools
  filter(richness > 0) |>
  #total nr of pools per cluster, across all sampling events
  mutate(p = length(unique(poolname)), .by = cluster) |>
  #number of pools per richness level in a given cluster at a given sampling event in a given year
  summarize(n = length(unique(poolname)), .by = c("cluster", "p", "sample", "year", "richness")) |> #, "author"
  #fraction of pools with given richness level
  mutate(fraction = n / sum(n), .by = c("cluster", "p", "sample", "year")) |> #, "author"
  #add desiccation data
  left_join(desiccation, by = join_by(cluster, year)) |>
  #add environmental heterogeneity
  left_join(heterogeneity, by= join_by(cluster)) 

#plot of fraction vs. nr of patches p or desiccation
test |>
  filter(sample == 1) |>
  ggplot() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x=desiccation, y=fraction, 
      col=as.factor(richness)) + 
  geom_point() +
  #facet_grid(as.factor(richness)~sample, scales = "free", labeller = label_both) +
  geom_smooth(lwd=0.5, method = lm, se=F, show.legend = F) + 
  labs(y = "fraction", col="richness")

ggplot(test) + 
  aes(x=p, y=desiccation) + 
  geom_point()

# ggsave("../figures/desiccation.pdf", width=4, height=3)
# ggsave("../figures/p.pdf", width=6, height=4)
# ggsave("../figures/desiccation-poly.pdf", width=6, height=4)

# Hard to see so check out the slopes of fraction vs. predictor (p or desiccation)
slopes <- test %>%
  nest_by(sample, richness, year) %>%
  expand_grid(predictor = c("p", "desiccation")) %>%
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

#Dominance of a given species?
dom <- counts |>
  mutate(richness = magna + longispina + pulex, 
         only_magna = (magna)*(1-longispina)*(1-pulex),
         only_longi = (1-magna)*(longispina)*(1-pulex),
         only_pulex = (1-magna)*(1-longispina)*(pulex)) |>
  select(!magna, !longispina, !pulex) |>
  pivot_longer(starts_with("only"), 
               names_to = "species", values_to = "value") |>
  #remove empty pools
  filter(richness > 0) |>
  #total nr of pools per cluster, across all sampling events
  mutate(p = length(unique(poolname)), .by = cluster) |>
  #fraction of pools with only species x
  summarize(fraction = sum(value) / length(unique(poolname)), 
            .by = c("cluster", "p", "sample", "species", "year")) |> #, "author"
  #add desiccation data
  left_join(desiccation, by = join_by(cluster))  

ggplot(dom) + 
  theme_bw() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x=desiccation, y=fraction, col=as.factor(year)) + 
  geom_point(show.legend = F) + 
  facet_grid(sample~species) +
  geom_smooth(lwd=0.5, method = lm, se=F, show.legend = F)
# ggsave("../figures/dom.pdf", width=6, height=4)


# Result 1: effect of mean desiccation: more single sp patches, fewer pairs, a bit more triplets
# Result 2: effect of nr of patches: fewer single sp patches 
# Interpretation: 
  # Homogeneity of environmental conditions? Or pH, conductivity, plants not most important? Heterogeneity similar across all clusters
  # Regional dominant? Not supported by the model? Pulex (and sometimes longi) according to the data?


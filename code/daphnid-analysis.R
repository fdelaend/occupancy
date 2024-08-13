
# READ DATA ----
island_measures <- read_csv("../data/Ebert/Frederik_data1/Island_measures_vers4.csv") |>
  select(island, group_100, group_200) |>
  rename(cluster = group_200) |>
  select(island, cluster)

counts          <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv") |>
  left_join(island_measures, by="island")

dessication     <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
  separate_wider_delim(poolname, delim = "-", names = c("island", "pool")) |>
  left_join(island_measures, by="island") |>
  summarise(dessication = mean(meandesi, na.rm = T), 
            p=n(), .by = "cluster") |>
  remove_missing()

plant_cover <- read_csv("../data/Ebert/Frederik_data1/plantcover_2013_2017.csv") |>
  remove_missing() |>
  # normalise across whole data set
  summarize(plants = mean(plants_rank), 
            .by = c("year", "sample", "island", "pool")) |>
  mutate(plants = (plants-mean(plants))/sd(plants))

pH_conductivity <- read_csv("../data/Ebert/Frederik_data1/MetapopData_pH_conductivity_1998_2017.csv") |>
  remove_missing() |>
  # compute mean data for when assessed multiple times
  summarize(pH = mean(pH), 
            log_conduct = mean(log10(conduct_uS)), 
            .by = c("year", "sample", "island", "pool")) |>
  # normalise across whole data set
  mutate(pH = (pH-mean(pH))/sd(pH),
         log_conduct = (log_conduct-mean(log_conduct))/sd(log_conduct)) 

ggplot(pH_conductivity) + 
  aes(x=year, y=log_conduct, col=pool) + 
  geom_point(show.legend = F) + 
  geom_smooth(method = "lm", show.legend = F, se=F, lwd=0.5)

# ANALYSE DATA ----
# compute heterogenity
heterogeneity <- pH_conductivity |>
  left_join(plant_cover, by = join_by(year, sample, island, pool)) |>
  #compute slope of environmental variable vs. time
  #first pivot so we don't need to do the same thing 3x
  pivot_longer(cols=c("pH", "log_conduct", "plants"), 
               values_to = "value", names_to = "variable") |>
  #now get the slope
  mutate(slope = cov(year,value)/var(year), 
         .by = c("variable", "sample", "island", "pool")) |>
  #add cluster ID
  left_join(island_measures, by = join_by(island)) |>
  #now get mean slope per cluster
  summarise(mean_slope = mean(slope, na.rm=T), 
            .by = c("sample", "variable", "cluster"))

  # environmental heterogeneity across pools within a cluster
  summarize(heterogeneity = sqrt(1/n()*sum((pH-mean(pH, na.rm=T))^2+
                                             (log_conduct-mean(log_conduct, na.rm=T))^2+
                                             (plants-mean(plants, na.rm=T))^2)), 
            .by = c("year", "sample", "cluster")) |>
  remove_missing()

test<-counts |>
  mutate(richness = magna + longispina + pulex) |>
  #remove empty pools
  filter(richness > 0) |>
  #number of pools per richness level in a given cluster at a given sampling event in a given year
  summarize(n = n(), .by = c("cluster", "richness", "sample", "year")) |> #, "author"
  mutate(fraction = n / sum(n), .by = c("cluster", "sample", "year")) |> #, "author"
  left_join(dessication, by = join_by(cluster)) |>
  #add environmental heterogeneity
  left_join(heterogeneity, by= join_by(year, sample, cluster)) 

test |>
  ggplot() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=p, y=fraction, col=heterogeneity) + 
  geom_point() + 
  facet_grid(richness~sample, scales = "free", labeller = label_both) +
  geom_smooth(aes(color=heterogeneity), lwd=0.5, 
              method = lm, se=F, show.legend = F)

ggsave("../figures/dessication.pdf", width=6, height=6)

# SUPPLEMENTS ----

#check out location of pools
pools           <- read_csv("../data/Ebert/Frederik_data1/Pools_coordinates_2017_vers7.csv") |>
  filter(core == 1) |>
  left_join(island_measures, by="island") |>
  select(name, island, pool, latitude_corr, longitude_corr, cluster) |>
  mutate(meanLat = mean(latitude_corr), #compute centroid
         meanLong = mean(longitude_corr), .by = c("cluster"))

ggplot(pools) +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x=latitude_corr, y=longitude_corr, col=as.factor(cluster)) +
  geom_point() +
  geom_point(aes(x=meanLat, y=meanLong), col="black")



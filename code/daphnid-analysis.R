
#read data ----
island_measures <- read_csv("../data/Ebert/Frederik_data1/Island_measures_vers4.csv") |>
  select(island, group_100, group_200) |>
  rename(cluster = group_200) |>
  select(island, cluster)

counts          <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv") |>
  left_join(island_measures, by="island")

pools           <- read_csv("../data/Ebert/Frederik_data1/Pools_coordinates_2017_vers7.csv") |>
  filter(core == 1) |>
  left_join(island_measures, by="island") |>
  select(name, island, pool, latitude_corr, longitude_corr, cluster) |>
  mutate(meanLat = mean(latitude_corr), #compute centroid
         meanLong = mean(longitude_corr), .by = c("cluster"))

dessication     <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
  separate_wider_delim(poolname, delim = "-", names = c("island", "pool")) |>
  left_join(island_measures, by="island") |>
  summarise(dessication = sum(meandesi, na.rm = T), 
            p=n(), .by = "cluster") |>
  remove_missing()

#check out location of pools
ggplot(pools) +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x=latitude_corr, y=longitude_corr, col=as.factor(cluster)) +
  geom_point() +
  geom_point(aes(x=meanLat, y=meanLong), col="black")

#analyse data ----
test<-counts |>
  mutate(richness = magna + longispina + pulex) |>
  #remove empty pools
  filter(richness > 0) |>
  #number of pools per richness level in a given cluster at a given sampling event in a given year
  summarize(n = n(), .by = c("cluster", "richness", "sample", "year")) |> #, "author"
  #year
  mutate(fraction = n / sum(n), .by = c("cluster", "sample", "year")) |> #, "author"
  left_join(dessication, by = "cluster") |>
  remove_missing()

test |>
  ggplot() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  aes(x=dessication, y=fraction, col=as.factor(year)) + 
  geom_point(show.legend = F) + 
  facet_grid(richness~sample, scales = "free", labeller = label_both) +
  geom_smooth(aes(color=as.factor(year)), lwd=0.5, 
              method = lm, se=F, show.legend = F)

ggsave("../figures/dessication.pdf", width=6, height=6)

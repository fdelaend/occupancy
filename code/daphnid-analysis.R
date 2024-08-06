
#read data ----
counts   <- read_csv("../data/Ebert/Frederik_data1/Daphnia_dynamics_1982_2017_2.csv")
pools    <- read_csv("../data/Ebert/Frederik_data1/Pools_coordinates_2017_vers7.csv") |>
  select(all_of(c("name", "island", "pool", "latitude_corr", "longitude_corr"))) |>
  mutate(meanLat = mean(latitude_corr), #compute centroid
         meanLong = mean(longitude_corr), .by = c("island"))
dessication <- read_csv("../data/Ebert/Frederik_data1/Data_hydroperiod_Means.csv") |>
  separate_wider_delim(poolname, delim = "-", names = c("island", "pool")) |>
  summarise(dessication = sum(mediandesi, na.rm = T), .by = "island")

#check out location of pools
ggplot(pools) +
  aes(x = latitude_corr, y = longitude_corr) +
  geom_point() +
  facet_wrap(vars(island), scales = "free") + 
  geom_point(aes(x = meanLat, y = meanLong), col="red")

# compute dispersion of pools and add number of pools to see if related
islandsDispersion <- pools |>
  summarise(dispersion = sqrt(mean((latitude_corr-meanLat)^2 + (longitude_corr-meanLong)^2)),
            nb = n(),
            .by = c("island"))
# yes, related, so probably not the best proxies of dispersal
ggplot(islandsDispersion) + 
  aes(x = nb, y = dispersion) + 
  geom_point()

#analyse data ----
test<-counts |>
  mutate(richness = magna + longispina + pulex) |>
  #remove empty pools
  filter(richness > 0) |>
  #number of pools per richness level on a given island at a given sampling event in a given year
  summarize(n = n(), .by = c("island", "richness", "sample", "year")) |> #, "author"
  #year
  mutate(fraction = n / sum(n), .by = c("island", "sample", "year")) |> #, "author"
  left_join(dessication, by = "island") |>
  remove_missing()

ggplot(test) +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=richness, y=fraction, col=dessication, 
      group=dessication) + 
  geom_line() +
  #geom_smooth(aes(group=dessication), se=F) +
  facet_grid(.~sample)

model <- lm(fraction ~ dessication + richness + dessication*richness,
            data=test|>filter(sample==1))

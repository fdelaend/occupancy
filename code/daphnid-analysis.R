source("tools-other.R") # function to fit glmms for priority effects

# Read and organise data ------------

## Cluster data: List of islands and clusters they belong to -------
cluster_measures <- read_csv("../daphnid-data/island-measures.csv") |>
  mutate(cluster = factor(group_100)) |>
  select(cluster, island) 

## Spatial location of pools and nr of surrounding pools within cluster --------- 
pool_coords <- read_csv("../daphnid-data/pool-coords.csv") |>
  # Pool coordinates 
  filter(!grepl("A", pool), !grepl("A", name)) |> #Only 1 observation for this island for the preemptive comp. analysis
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

## Plot of pool locations ------
pool_map <- pool_coords %>%
  st_as_sf(coords = c("longitude_corr", "latitude_corr"), crs = 4326) %>%
  mutate(cluster = factor(cluster))

# Finland outline
finland <- ne_countries(country = "Finland", scale = "medium", returnclass = "sf")

# Use a projected CRS suitable for Finland
pool_map_proj <- st_transform(pool_map, 3067)
finland_proj <- st_transform(finland, 3067)

# Main map extent: buffer around your points
bbox <- st_bbox(
  st_buffer(st_union(pool_map_proj), dist = 200)   # buffer in meters
)

main_map <- ggplot() +
  geom_sf(data = finland_proj, fill = "grey95", color = "grey70") +
  geom_sf(data = pool_map_proj, aes(color = cluster), size = 1, alpha = 0.9) +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE) +
  scale_colour_manual(values = cbPalette) + 
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  labs(color = "Cluster", x = NULL, y = NULL)

# Inset map of Finland with your region marked
region_box <- st_as_sfc(st_bbox(pool_map_proj)) %>% st_buffer(10000)

inset_map <- ggplot() +
  geom_sf(data = finland_proj, fill = "grey95", color = "grey70", linewidth = 0.3) +
  geom_sf(data = region_box, fill = NA, color = "black", linewidth = 1) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       color = "black",
                                       linewidth = 0.8),
        plot.margin = margin(4, 4, 4, 4))

# Combine main map + inset
ggdraw(main_map) +
  draw_plot(inset_map, 
            x = 0.5, y = 0.15, width = 0.28, height = 0.2)

ggsave(filename = "../figures/map.pdf", 
       width=5, height = 6, device = "pdf")

## Daphnid data -------
# Read Daphnid data and join cluster info and nr of pools within "b"
daphnids          <- read_csv("../daphnid-data/daphnid-counts.csv") |>
  filter(!grepl("A", pool), !grepl("A", poolname)) |> #Only 1 observation for this island for the preemptive comp. analysis
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

### Plot for talk ------
# Select the relevant columns
species_df <- data_priority %>%
  select(
    magna_spring, magna_summer_alone,
    longispina_spring, longispina_summer_alone,
    pulex_spring, pulex_summer_alone
  )

# Convert to long format
long_df <- bind_rows(
  tibble(
    species = "magna",
    spring = data_priority$magna_spring,
    summer_alone = data_priority$magna_summer_alone
  ),
  tibble(
    species = "longispina",
    spring = data_priority$longispina_spring,
    summer_alone = data_priority$longispina_summer_alone
  ),
  tibble(
    species = "pulex",
    spring = data_priority$pulex_spring,
    summer_alone = data_priority$pulex_summer_alone
  )
)

# Calculate P(summer_alone = 1 | spring)
prob_df <- long_df %>%
  group_by(species, spring) %>%
  summarise(
    prob = mean(summer_alone == 1, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Plot
ggplot(prob_df,
       aes(x = factor(spring),
           y = prob)) +
  geom_col(width = 0.7) +
  facet_wrap(~ species) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    x = "Spring presence (0/1)",
    y = "P(Summer alone = 1)",
    title = "Probability of summer-alone occurrence conditional on spring occurrence"
  ) +
  theme_classic()

ggsave(filename = "../figures/priority-talk.pdf", 
       width=5, height = 3, device = "pdf")

### Model fitting --------
all_species <- c("magna", "longispina", "pulex")

models_priority <- tibble(focal = all_species) |>
  mutate(model = map(focal, 
                     ~ fit_priority(focal = .x, adapt_delta = 0.999,
                                      data = data_priority))) 

saveRDS(models_priority, "../data/models_priority.RDS")
models_priority <- readRDS("../data/models_priority.RDS")

### Model plotting --------
# get out the typical names of the parameters we want to plot
colnames(as.matrix(models_priority$model[[1]]))
# do the plotting: Fig 4
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
  scale_fill_manual(values = cbPalette) + 
  aes(x = value, y = focal, fill = parameter) + 
  stat_halfeye(alpha = 0.6, slab_color = NA) +
  #facet_wrap(~ parameter, scales = "free_x") +
  labs(x = "parameter value", y= "focal species") +
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

### Model printing -----
print_model(models_priority, focal)
                               
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
# do the plotting: Fig 5
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
  scale_fill_manual(values = cbPalette) + 
  aes(x = value, y = as.factor(b), fill = parameter) + 
  stat_halfeye(alpha = 0.6, slab_color = NA) +
  facet_grid(sample~parameter, 
             scales = "free") +
  theme_bw() +
  labs(y = "distance (m)", x = "parameter value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(filename = "../figures/proportion.pdf", 
       width=12, height = 3, device = "pdf")

#### Fixed effect: Fig 6 -------
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
       width=6, height = 3, device = "pdf")

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

### Model printing -------------
print_model(models_prop, sample, b)


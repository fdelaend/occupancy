
## load the data
counts_env <- readRDS("../data/dataDaphnids.rds") # To change
# variables are: 
# year: year of sampling
# sample: spring (1), or summer sample (2)
# island: unique island ID
# pool: unique pool ID
# cluster: unique cluster ID. A cluster is a collection of pools lying within a given distance from each other
# richness: Daphnia species richness
# pools_dens: number of pools per running m of a cluster's perimeter
# desiccation_dynamic: estimates of pool desiccation rate
# pH: measured pH

# Check correlations among potential predictors of richness
counts_env |> 
  select(year, sample, pools_dens, desiccation_dynamic, pH) |>
  cor(use = "pairwise.complete.obs")
# desiccation_dynamic and year seem to be correlated and 
# therefore I do not put them in the same model. 
# It is well known that pools dry more often in recent years,
# and so I assume that a year effect would mostly be a dryness effect.
# I further drop pH from the list of predictors because we don't have 
# pH measurements for 90% of the data. 
# Finally, I focus on the summer sample because most representative. 

# Fit the models 
# Mixed GLM
model1 <- glmer(richness ~ pools_dens + desiccation_dynamic + (1|cluster),
                data = counts_env |> filter(sample == 2), 
                family = poisson(link = "log"))
summary(model1)
# Desiccation has positive effect on richness. No effect of pool density
# But maybe driven by ponds at richness = 3?
ggplot(counts_env |> filter(sample == 2)) + 
  aes(y = as.factor(richness), x = desiccation_dynamic) + 
  geom_boxplot()

# Regular GLM
model2 <- glm(richness ~ pools_dens + desiccation_dynamic, 
              data = counts_env |> filter(sample == 2), 
              family = poisson(link = "log"))
summary(model2)
# Results are different: now pools density has positive effect, 
# while desiccation has negative effect. Trends look dodgy though, 
# and seem to be driven by data at richness = 2
ggplot(counts_env |> filter(sample == 2)) + 
  aes(y = as.factor(richness), x = pools_dens) + 
  geom_boxplot()

# Problem might be that there are only 31 observations with 
# richness = 3 (out of the 30000). 
# What happens if we remove these ponds? 
# Mixed GLM
model1 <- glmer(richness ~ pools_dens + desiccation_dynamic + (1|cluster),
                data = counts_env |> filter(sample == 2, richness < 3), 
                family = poisson(link = "log"))
summary(model1)
# Same result as before

# Regular GLM
model2 <- glm(richness ~ pools_dens + desiccation_dynamic, 
              data = counts_env |> filter(sample == 2, richness < 3), 
              family = poisson(link = "log"))
summary(model2)
# Same result as before

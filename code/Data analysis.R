
#read data ----
individuals   <- read_csv("../data/doi_10_5061_dryad_8gp03__v20160317/Sample_information.csv")
groomingViola <- read_csv("../data/doi_10_5061_dryad_8gp03__v20160317/Viola_1yrgroom.csv")
groomingMica  <- read_csv("../data/doi_10_5061_dryad_8gp03__v20160317/Mica_1yrgroom.csv")
dViola        <- sum(groomingViola[,c(2:ncol(groomingViola))])
dMica         <- sum(groomingMica[,c(2:ncol(groomingMica))])
dietViola     <- read_csv("../data/doi_10_5061_dryad_8gp03__v20160317/Viola_1yr_diet_BCdistance.csv")
dietMica      <- read_csv("../data/doi_10_5061_dryad_8gp03__v20160317/Mica_1yr_diet_BCdistance.csv")
hist(unlist(dietViola[,-1]))
hist(unlist(dietMica[,-1]))

abundances    <- read_csv("../data/doi_10_5061_dryad_8gp03__v20160317/Abundances_species_MetaPhlAn2.csv") %>%
  select(-Species) 
abundances <- abundances %>%
  filter(rowSums(abundances)>0.1) #kick out species that occur almost nowhere
richness <- colSums(abundances>0)

#analyse data ----
# nr of patches and intensity of interactions
baboons <- individuals %>%
  mutate(richness = richness[base::match(individuals$`Subject ID`, names(richness))]) %>%
  drop_na(richness)

ggplot(baboons) + 
  theme_bw() +
  aes(x=richness, col=`Social group`) + 
  geom_density()

ggplot(abundances %>% pivot_longer(cols=everything(), names_to="baboon",
                                   values_to = "abundance")) + 
  theme_bw() +
  aes(x=log10(abundance), col=baboon) + 
  geom_density(show.legend=F)

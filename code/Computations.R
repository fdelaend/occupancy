
# CALCULATIONS ------
data <- expand_grid(n = c(3), meanA = c(0, 0.4, 0.8), #, 6
                    d = seq(-8,-4, length.out=3), #6
                    sdA = 0, p = 100, rep = c(1:3)) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters
  mutate(d = 10^d) %>%
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n))) %>% 
  mutate(A = map(A, ~make_symmetric(.x))) %>% #make symmetric and ditch diagonal
  mutate(A = map(A, ~ set_diagonal(A=.x, d=1))) %>% #set self to 1
  mutate(feasibility = map_dbl(A, ~feasibility(.x))) %>% #compute feasibility
  mutate(A = pmap(., make_block_diagonal)) %>% #make A spatial
  mutate(R = pmap(., make_R_spatial)) %>% #make spatial R (note that the emigration is subtracted later so this is the real local R)
  mutate(D = pmap(., make_D)) %>% #make dispersal matrix
  mutate(N0 = map(R, ~ .x*0+extinctionThreshold)) %>% #set initial conditions
  #Simulate the network
  mutate(NHat = pmap(., get_NHat)) %>%
  #Add R for later analysis
  mutate(NHat = map2(NHat, R, ~ .x %>% mutate(R=.y))) %>%
  #Analyse: 
  #0/cv of R within a patch
  mutate(cvR = map(NHat, ~ .x %>% group_by(location) %>% summarize(cvR = sd(R)/mean(R)))) %>%
  #1/nr of patches where a subset m<= n persists
  mutate(nrPatchesM = pmap(., function(NHat, n,...) {
    NHat %>% mutate(present = density>extinctionThreshold) %>%
      group_by(location) %>%
      summarize(m = as_factor(sum(present))) %>%
      mutate(m = fct_expand(m, as.character(c(1:n)))) %>%
      group_by(m) %>%
      count(m, .drop=F, name="nrPatches")})) %>%
  #1.1/nr of patches in which all n species persist
  mutate(nrPatchesN = map2_dbl(nrPatchesM, n, ~ (.x %>% filter(m==.y))$nrPatches)) %>%
  mutate(propPatchesN = nrPatchesN/p) %>%
  #2/total density across all patches and species
  mutate(NTot = map_dbl(NHat, ~ sum(.x["density"]))) %>%#total density across all patches
  #3/mean density of each species across all patches  
  mutate(meanAcross = pmap(., function(NHat,...) {
    NHat %>% 
      group_by(sp) %>%
      summarize(meanNAcross = mean(density))})) %>%
  #4/mean density at each site, across all species
  mutate(meanWithin = pmap(., function(NHat,...) {
    NHat %>% 
      group_by(location) %>%
      summarize(meanNWithin = mean(density))})) %>%
  #5/total of local interactions experienced by every species: R - AN, 
  # after having set diag(A) to zero (because we only want interactions with heterospecifics)
  mutate(ANoti = map(A, ~ set_diagonal(A=.x, d=0))) %>% #set self to 0
  mutate(X = map2(NHat, ANoti, ~ c(.x$R)-.y%*%c(.x$density))) %>%
  # add info to NHat solution
  mutate(NHat = map2(NHat, X, ~ .x %>% mutate(X=.y[,1]))) %>%
  mutate(XMean = map(NHat, ~ .x %>% group_by(location) %>% 
                       summarize(XMean=mean(X^2))))

# PLOTS ------
## main plot ----
ggplot(data) + 
  scale_color_gradient(low = "yellow", high = "red") +
  scale_linetype_discrete(rep("solid", 100)) +
  theme_bw() +
  aes(x=log10(d), y=propPatchesN, col=meanA) + 
  geom_point() +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy (fraction)", col="a") +
  geom_line(aes(x=log10(d), y=feasibility, col=meanA,
                group=interaction(meanA, rep))) + 
  facet_grid(cols=vars(n))

ggsave(paste0("../figures/feas.pdf"), width=6, height = 2, 
       device = "pdf")

## total biomass ----
ggplot(data %>% 
         mutate(NTotPred = p*get_N_total(mean_a=meanA, d=1, l=n, r=1))) + 
  scale_linetype_discrete(rep("solid", 100)) +
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=log10(d), y=NTot, col=meanA) + 
  labs(x=expression(paste("log"[10],"(d)")), 
       y="total biomass \n across network", col="a") +
  geom_point() + 
  facet_grid(cols=vars(n))

ggsave(paste0("../figures/biomass.pdf"), width=5, height = 2, 
       device = "pdf")

## the number of patches occupying m<=n species in case there is no dispersal:
dataNoDisp <- data %>%
  select(n, A, meanA, d, p, rep, nrPatchesM) %>%
  filter(d==min(d)) %>%
  unnest(nrPatchesM) %>%
  mutate(m=as.numeric(as.character(m))) %>%
  mutate(fractionPatches = nrPatches/p) %>%
  mutate(fractionPatchesPredictedOne = 
           map2_dbl(A, n, ~(1-feasibility(.x[c(1:2),c(1:2)]))^(.y-1))) %>%
  filter(m==1)
  
ggplot(dataNoDisp) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=m, y=fractionPatchesPredicted, 
      col=meanA) + 
  geom_point() + 
  facet_grid(cols=vars(n)) #+ 
  #geom_abline(slope=1, intercept = 0)

## Does a well-mixed system behave as a single system?
# No
# But switching off dispersal gives same pattern,
# so probably only due to the fact that we have many patches, 
# some of which contain both, others only 1, still others only 2. 

  
## Can we predict the total density across the whole network 
## with a model that acts as if all were a single system with summed Rs as the R?
ggplot(data) + 
  scale_color_gradient(low = "yellow", high = "red") +
  theme_bw() +
  aes(x=NsingleTotal, y=NTot, col=d) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
#No!

# LEFTOVERS ---------

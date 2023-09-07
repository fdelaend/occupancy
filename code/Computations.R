extinctionThreshold <- 1e-3

# CALCULATIONS ------
data <- expand_grid(n = c(3, 6), meanA = c(0, 0.2, 0.4, 0.8), 
                    d = seq(-6,-4, length.out=6),
                    sdA = 0, p = 50, rep = c(1:3)) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters
  mutate(d = 10^d) %>%
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n))) %>% 
  mutate(A = map(A, ~make_symmetric(.x))) %>% #make symmetric and ditch diagonal
  mutate(A = map(A, ~ set_diagonal(A=.x, d=1))) %>% #set self to 1
  #compute feasibility
  mutate(feasibility = map_dbl(A, ~feasibility(.x))) %>%
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
  #1/nr of feasible patches (aka where all present)
  mutate(nrFeas = pmap_dbl(., function(NHat, n,...) {
    NHat %>% mutate(present = density>extinctionThreshold) %>%
      group_by(location) %>%
      summarize(nrSp = sum(present)) %>%
      filter(nrSp==n) %>% nrow()})) %>%
  mutate(propFeas = nrFeas/p) %>%
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
ggplot(data) + 
  scale_color_gradient(low = "yellow", high = "red") +
  scale_linetype_discrete(rep("solid", 100))+
  theme_bw() +
  aes(x=log10(d), y=propFeas, col=meanA) + 
  geom_point() +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy (fraction)", col="a") +
  geom_line(aes(x=log10(d), y=feasibility, col=meanA,
                group=interaction(meanA, rep))) + 
  facet_grid(cols=vars(n))

ggsave(paste0("../figures/feas.pdf"), width=4, height = 3, 
       device = "pdf")

## Does a well-mixed system behave as a single system?
# No
# But switching off dispersal gives same pattern,
# so probably only due to the fact that we have many patches, 
# some of which contain both, others only 1, still others only 2. 
ggplot(data %>% filter(d==10^-6)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  theme_bw() +
  aes(x=`1`, y=`1single`, col=meanA) + 
  geom_point() + 
  labs(x="total density sp 1 across network", 
       y="density of sp 1 in one mean system") + 
  geom_abline(slope = 1, intercept = 0)

## Do patches with two species contain more density than pathces with one?
density_per_patch <- data %>%
  select(all_of(c("n", "meanA", "d", "totalPerSite"))) %>%
  unnest(totalPerSite)

ggplot(density_per_patch) + 
  scale_color_gradient(low = "yellow", high = "red") +
  theme_bw() +
  aes(x=nrSp, y=NTotal_loc, col=meanA) +
  geom_point()
  
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
  #5/model equilibrium as if all were a single system with summed Rs as the R
  mutate(Rsingle = pmap(., function(R, n, p,...) {
    tibble(R) %>% 
      mutate(sp = rep(c(1:n), p)) %>%
      group_by(sp) %>%
      summarize(R = sum(R)) %>% #total across whole network
      select(R) %>%
      as_vector()})) %>%
  mutate(Nsingle = map2(meanA, Rsingle, ~ solve(diag(length(.y))+.x-diag(length(.y))*.x)%*%.y)) %>%
  #get total density as such obtained
  mutate(NsingleTotal = map_dbl(Nsingle, ~ ifelse(prod(.x>0), sum(.x), NA))) %>%
  #just pick one species 
  mutate(`1single` = map_dbl(Nsingle, ~ifelse(prod(.x>0), .x[1], NA)))
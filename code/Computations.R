
# CALCULATIONS ------
data <- expand_grid(n = c(6), meanA = c(0.2, 0.8), #, 6
                    d = seq(-8,-4, length.out=2), #6
                    sdA = 0, p = 100, rep = c(1:2)) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
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
  #1/nr of patches where a subset m<= n persists
  mutate(nrPatchesM = pmap(., function(NHat, n,...) {
    NHat %>% mutate(present = density>extinctionThreshold) %>%
      group_by(location) %>%
      summarize(m = as_factor(sum(present))) %>%
      mutate(m = fct_expand(m, as.character(c(1:n)))) %>%
      group_by(m) %>%
      count(m, .drop=F, name="nrPatches")})) %>%
  #2/nr of patches in which all n species persist
  mutate(nrPatchesN = map2_dbl(nrPatchesM, n, ~ (.x %>% filter(m==.y))$nrPatches)) %>%
  mutate(propPatchesN = nrPatchesN/p) %>% #proportion
  #3/total density across all patches of a species
  mutate(NTotalK = pmap(., function(NHat,...) {
    NHat %>% 
      group_by(sp) %>%
      summarize(NTotalK = sum(density))}))

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

## showcase accuracy of f(m) -----
dataNoDisp <- data %>%
  select(n, A, meanA, d, p, rep, nrPatchesM) %>%
  filter(d==min(d), meanA>0) %>%
  unnest(nrPatchesM) %>%
  mutate(m=as.numeric(as.character(m))) %>%
  mutate(fractionPatches = nrPatches/p) %>%
  mutate(fractionPatchesTest = pmap_dbl(., get_fraction_m))

ggplot(dataNoDisp) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=fractionPatches, y=fractionPatchesTest, col=meanA) + 
  geom_point() + 
  geom_abline(slope=1, intercept = 0)

ggplot(dataNoDisp %>% mutate(meanA2 = meanA) %>%
         unite("it", rep, meanA2)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  geom_point(aes(x=m, y=fractionPatches, col=meanA), alpha=0.5) + 
  aes(x=m, y=fractionPatchesTest, col=meanA, lty=as_factor(it)) + #,,  
  geom_line(show.legend = F) + 
  facet_grid(cols=vars(n)) #+ 
  #geom_abline(slope=1, intercept = 0)

## Analytical predictions -----
#Case where i is extinct w/o dispersal
Prob1 <- expand_grid(meanA = c(0.01, 0.02, 0.6, 0.8), n=c(4, 6), sampleSize = 100) %>%
  mutate(NTotalValues = map(n, ~get_N_total(n=c(1:(.x)))), #total across sp; a RV
         mProbs = map2(n, meanA, ~make_distribution(n=.x, meanA=.y)),
         NMeanK = map2_dbl(NTotalValues, mProbs, ~sum(.x*.y))) %>% #mean across sites (intermediate result for NTotalK)
  expand_grid(d=seq(-6,-4, length.out=6), p=100) %>%
  mutate(NTotalK = NMeanK*p/n) %>%
  mutate(d=10^d) %>%
  mutate(NTotal = pmap(., function(NTotalValues, sampleSize, mProbs, ...) { #sample NTotal
    sample(x=NTotalValues, size=sampleSize, prob = mProbs, replace=T)})) %>%
  mutate(ri = map(sampleSize, ~sample(seq(1e-3,4,length.out=100), size=.x, 
                                      prob=dnorm(seq(0,4,length.out=100), 1, 0.7),
                                      replace = T))) %>%#sample from truncated normal
  mutate(r = map(ri, ~mean(.x)),  #compute mean and mean of inverse
         rInv = map(ri, ~mean(1/.x))) %>%
  mutate(m= map2(mProbs, sampleSize, ~sample(x=c(1:length(.x)), size=.y, prob=.x, replace=T))) %>% #sample m's
  mutate(NiStrong = pmap(., function(d, NTotalK, meanA, NTotal, ri, ...){
    d*NTotalK/(meanA*NTotal-ri)})) %>% #compute Ni for strong interactions
  mutate(NiWeak = pmap(., get_density_weak)) %>% #compute Ni for weak interactions
  mutate(Ni = pmap(., function(meanA, NiWeak, NiStrong,...){NiWeak*(meanA<0.5)+NiStrong*(meanA>=0.5)})) %>%#get proper Ni according to meanA
  mutate(ProbNi = pmap_dbl(., function(n, Ni, sampleSize, ...){ #prob 
                           (sum(Ni>extinctionThreshold)/sampleSize)^n})) %>%
  mutate(ProbNi0Zero = pmap_dbl(., function(n, mProbs, ...){
    sum(mProbs[1:(n-1)]*(1-c(1:(n-1))/n)) })) 

ggplot(Prob1) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=log10(d), y=ProbNi, col=meanA) + 
  geom_point() + 
  facet_grid(cols=vars(n)) #+ 

# LEFTOVERS ---------


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

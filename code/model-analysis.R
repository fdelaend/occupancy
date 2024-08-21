
# SIMULATIONS ----
## Simulations ------
Sims <- expand_grid(n = c(6), meanA = c(0.2, 0.4, 0.8, 1.2), #, 6; 0.4, 0.4, 0.6, 
                    d = seq(-6,-4, length.out=6), vary=c(0, 0.1), k=c(1, 1.5),
                    cvA = c(0, 0.5), p = 100, rep = c(1:5)) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters
  mutate(d = 10^d,
         sdA = cvA * meanA) %>%
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% #make symmetric and ditch diagonal; set self to 1
      make_symmetric() %>% set_diagonal(d=1))) %>% 
  mutate(A = pmap(., make_block_diagonal)) %>% #make A spatial
  mutate(R = pmap(., make_R_spatial)) %>% #make spatial R (note that the emigration is subtracted later so this is the real local R)
  mutate(regularD = pmap(., make_D)) %>% #make classic dispersal matrix
  mutate(coords = map(p, ~randomCoords(nPatch=.x)), #generate coordinates for every patch
         distances = map(n, ~seq(4,0.5, length.out=.x))) %>% #generate characteristic distances of the species
  mutate(exponentialD = map2(coords, distances, ~dispMatrixCommunityExp(coords=.x, #Make a D matrix for exponential decay
                                                              dispDistanceVector=.y)),
         exponentialD = map2(exponentialD, d, ~rescale_D(D=.x, d=.y))) %>%
  mutate(N0 = map(R, ~ .x*0+extinctionThreshold)) %>% #set initial conditions
  pivot_longer(cols = c("regularD","exponentialD"), names_to = "dispType", values_to = "D") %>% # pivot the two sorts of D matrices to make them a factor
  #Simulate the network
  mutate(NHat = pmap(., get_NHat)) %>%
  #1/ Summarize: per m, compute the nr of patches and total biomass of an average patch
  mutate(summaryM = pmap(., function(NHat, R, n,...) {
    NHat %>% mutate(present = density>extinctionThreshold) %>%
      mutate(R=R) %>%
      group_by(location) %>%
      summarize(m = as_factor(sum(present)),
                NTotal = sum(density),
                meanRPer = sum(R*present)/sum(present)) %>%
      mutate(m = fct_expand(m, as.character(c(1:n)))) %>%
      group_by(m, .drop=F) %>%
      summarize(nrPatches = n(),#nr of patches with m sp.
                NTotal = mean(NTotal),#total biomass in a patch with m sp.
                meanRPer = mean(meanRPer)) %>% #mean r of persisting sp.
      ungroup()})) %>%
  #2/proportion of patches in which all n species persist
  mutate(propPatchesN = 1/p*map2_dbl(summaryM, n, ~ (.x %>% filter(m==.y))$nrPatches)) %>%
  #3/total density across all patches of a species
  mutate(NTotalK = pmap(., function(NHat,...) {
    NHat %>% 
      group_by(sp) %>% #sum across patches for every sp.
      summarize(NTotalK = sum(density))})) #so you can check it's comparable across sp

#saveRDS(Sims, "../data/data.rds")
#Sims <- readRDS("../data/data.rds")

## Recover the distribution of the growth rates --------
Rs     <- Sims %>% filter(d==min(d), meanA==0.8, rep==1, k==1) %>% select(n, R, NHat) %>% 
  mutate(NHat = map2(R, NHat, ~ .y %>% mutate(R=.x))) %>%
  #filter(n==6) %>%
  select(-R) %>%
  unnest(NHat) %>%
  select(-density, -sp)
pdfRs  <- density(Rs$R, from=0)
meanR  <- sum(pdfRs$x*pdfRs$y)/sum(pdfRs$y)#grant mean of R

# INTERMEDIATE PREDICTIONS in absence of dispersal ----
## Make the predictions of f(m) and Ntotal -----
IntPred <- Sims %>% 
  # Select Sims w/o dispersal (select regularD although the same for both dispersal types),  
  # diffuse competition, and 
  # regional equivalence, and make predictions
  filter(d==min(d), dispType=="regularD", meanA<1, cvA==0, k==1, vary==0) %>%
  # only select useful bits
  select(d, n, meanA, p, rep, summaryM, NTotalK, vary, cvA, k) %>%
  mutate(NTotalK = map_dbl(NTotalK, ~mean(.x$NTotalK))) %>% #mean regional density across sp
  mutate(summaryM = pmap(., function(summaryM, meanA, n, ...) { #add predictions
    summaryM %>%
    mutate(m=as.numeric(as.character(m))) %>%
    rowwise() %>%
    mutate(fractionPatchesPredicted = get_fraction_m(meanA=meanA, m=m, n=n)) %>% #predicted total density in a patch of m persisting sp
    ungroup() %>%
    mutate(meanRPerPredicted = get_RMeanM(a=meanA, m=m, n=n, r=meanR),#predicted mean r of persisting sp
           meanRExcPredicted = (-meanRPerPredicted*m+n*meanR)/(n-m),#predicted mean r of excluded sp
           NTotalPredicted = get_N_total(meanA=meanA, n=m, r=meanRPerPredicted))})) %>% #predicted total density for a patch with m species
  mutate(NTotalKPredicted = p/n*map_dbl(summaryM, ~sum(.x$fractionPatchesPredicted * .x$NTotalPredicted))) #Total regional density for a species

## Unnest and only keep those variables that are for now relevant ---
SelIntPred <- IntPred %>%
  unnest(summaryM) %>%
  select(all_of(c("d", "meanA", "n", "p", "rep", "m", "fractionPatchesPredicted", 
                  "cvA", "vary", "k")))

## Add predictions of f(m) to simulations -----
IntPredsAndSims <- Sims %>%
  unnest(summaryM) %>%
  mutate(m=as.numeric(as.character(m))) %>%
  mutate(fractionPatches = nrPatches/p) %>% #get simulated f(m)
  select(all_of(c("d", "dispType", "meanA", "n", "p", "rep", "m", "fractionPatches", 
                  "cvA", "vary", "k"))) %>% #only keep relevant variables
  left_join(SelIntPred, by=c("d", "meanA", "n", "p", "rep", "m", "cvA", "vary", "k"),
            multiple = "all") %>% #join with predictions; 
  #makes sure that predictions are NA where dispersal is > min(d),
  #or where there is no regional equivalence
  group_by(d, dispType, meanA, n, p, m, cvA, vary, k) %>% 
  #compute summary stats of predicted and simulated f(m) 
  #(mean and sd)
  summarise(meanProb = mean(fractionPatches),
            sdProb = sd(fractionPatches),
            meanProbPred = mean(fractionPatchesPredicted),
            sdProbPred = sd(fractionPatchesPredicted))

## plot simulations vs. predictions of f(m); only take k==1 and vary==0 -----
IntPredsAndSims %>% 
  mutate(d=round(log10(d),1)) %>%
  filter(k==1, dispType=="regularD", d %in% c(-6, -5.2, -4.4), vary==0) %>% 
  ggplot() + 
  theme_bw() +
  scale_linetype_manual(values=rep("dashed", 1000)) + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_line(aes(x=m, y=meanProb, col=meanA, group=meanA), 
            alpha=0.4) + 
  geom_errorbar(aes(x=m, ymin=meanProb-sdProb, 
                    ymax=meanProb+sdProb, col=meanA, 
                    width = 0.1)) +
  geom_line(aes(x=m, y=meanProbPred, col=meanA, 
                lty=as_factor(meanA)), show.legend = F, lwd=1) + 
  facet_grid(cvA~d, 
             labeller = label_bquote(rows=paste("cv(", a[ij],")=", .(cvA)),
                                     cols=paste("log(D)= ", .(d)))) +
  labs(x="nr of persisting species, m", 
       y="fraction of patches with m species, f(m)", 
       col=expression(paste(bar(a))))

ggsave(paste0("../figures/fm.pdf"), width=4.5, height = 3, 
       device = "pdf")

## Plot predicted vs. simulated NtotalK ----
NtotalK <- ggplot(IntPred %>% group_by(meanA) %>% 
         summarise(meanNTotalK = mean(NTotalK),
                   sdNTotalK = sd(NTotalK),
                   meanNTotalKPredicted = mean(NTotalKPredicted),
                   sdNTotalKPredicted = sd(NTotalKPredicted))) + 
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=meanA, y=meanNTotalK) + 
  geom_point() + 
  geom_errorbar(aes(x=meanA, ymin=meanNTotalK-sdNTotalK, 
                    ymax=meanNTotalK+sdNTotalK, 
                    width = 0.05)) + 
  geom_line(aes(x=meanA, y=meanNTotalKPredicted)) +
  labs(x="interaction strength, a", 
       y=expression(paste("total density, ", Sigma[{k}],"N"[{"0,i"}]^{(k)}))) 
#labs(x=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"simulated")), 
#     y=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"analytical"))) +
#facet_grid(cols=vars(cvA), rows=vars(vary), labeller = label_both)

ggsave(paste0("../figures/Nk.pdf"), width=4, height = 3, 
       device = "pdf")

# PREDICTIONS OF PATCH OCCUPANCY -----
## Make the predictions for when all assumptions are met
Predictions <- IntPred %>%
  filter(rep==1) %>% #Select 1st replicate
  select(-c("rep", "k", "cvA", "vary", "d")) %>% #Remove irrelevant variables
  mutate(sampleSize = 1000) %>% #set sample size for probability calculations
  mutate(samples = pmap(., function(summaryM, sampleSize, meanA, NTotalKPredicted, p, ...){
    tibble(m = sample(x=summaryM$m, size=sampleSize, prob = summaryM$fractionPatchesPredicted, replace=T)) %>% #sample m according to f(m). Every sample is a hypothetical patch      
      left_join(summaryM, by="m", multiple="all") %>% #get variables that match the sampled m
      select(all_of(c("m", "NTotalPredicted", "meanRPerPredicted", "meanRExcPredicted"))) %>%
      rowwise() %>%       
      #sample 1 growth rate per hypo patch, from distribution of R such that IGR<0 (for riExc) or >0 (for riPer) 
      mutate(riExc=sample_ri(samplesize=1, PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="above"),
             riPer=sample_ri(samplesize=1, PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="below")) %>%
      ungroup() %>%
      mutate(N0i=get_N0i(a=meanA, n=m, r=meanRPerPredicted, ri=riPer), #predict density of persisting sp in absence of dispersal
             N1iExc=get_N1iExc(NTotalK=NTotalKPredicted, ri=riExc, a=meanA, NTotalPredicted), #density contribution per unit of dispersal, in case of exclusion w/o dispersal
             NegIGR=1/get_N1iExc(NTotalK=1, ri=riExc, a=meanA, NTotalPredicted), #negative IGR
             rhoi = NTotalKPredicted/p*(p-1)/N0i)})) %>% #ratio of total to local density
  mutate(means = map(samples, ~ .x %>% group_by(m) %>% 
                       summarize(meanN1iExc=mean(N1iExc, na.rm=T), #mean across patches of Ni1 in case of exclusion w/o dispersal. Because we sample 1 species per patch, this is the same as taking a mean across species
                                 rho = mean(rhoi, na.rm=T)))) %>% #same type of mean, but now of rhoi 
  mutate(samples = map2(samples, means, ~.x %>% #add means to samples
                          left_join(.y, by="m"))) %>%
  expand_grid(d = seq(-6,-4, length.out=60)) %>% #Create a gradient of dispersal values
  mutate(d=10^d) %>%
  mutate(samples = pmap(., function(samples, d, meanA, n, p, ...){ samples %>%
      mutate(N1iPer = get_N1iPer(a=meanA, n=n, m=m, rho=rho, #N1i when i persists w/o disp.
                                 rhoi=rhoi, meanN1iExc=meanN1iExc, p=p), 
             NiExc = d*N1iExc, #Ni when i is excluded w/o disp.
             NiPer = N0i+d*N1iPer)})) %>% #Ni when i persists w/o disp.
  mutate(probExc = map_dbl(samples, ~sum(.x$NiExc>extinctionThreshold)/length(.x$NiExc)), #Compute probabilities that > threshold
         probPer = map_dbl(samples, ~sum(.x$NiPer>extinctionThreshold)/length(.x$NiPer))) %>%
  mutate(summaryM = map2(summaryM, n, ~.x %>%
                             mutate(probN0iExt = fractionPatchesPredicted*(1-m/.y)))) %>%#proba that Ni is absent from patches with m sp (w/o disp)
  mutate(probN0iExt = map_dbl(summaryM, ~sum(.x$probN0iExt)), #overall proba across all patches
         probN0iPer = 1-probN0iExt,
         prob = (probExc*probN0iExt + probPer*probN0iPer)^n) %>% #grant prob
  select(-summaryM)

## Illustrate predictions of NtotalK and negative IGR ----
NegIGR <- Predictions %>% 
  filter(d==min(d)) %>%
  select(all_of(c("meanA", "samples"))) %>%
  unnest(samples) %>%
  ggplot() +
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=NegIGR, col=meanA, group=meanA) + 
  geom_density(show.legend = F) + 
  labs(x="negative invasion growth rate", 
       y="probability density", col="a")

## Plot predictions of prob, conditional on i persisting/excluded w/o disp. ----
probExc <- ggplot(Predictions) + 
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=log10(d), y=probExc, col=meanA, group=meanA) + 
  geom_line() +
  #geom_point() +
  labs(x=expression(paste("log"[10],"(d)")), 
       y=expression(paste("P(N"[i],">0 | N"[i,0],"<0)")),
       col="a")

case1 <- grid.arrange(NtotalK, NegIGR, probExc, 
                             ncol=3, widths=c(1,1,1.4)) 

ggsave(paste0("../figures/case1.pdf"), case1, 
       width=8, height = 3, device = "pdf")  

## Plot Case 2
probExc <- ggplot(Predictions) + 
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=log10(d), y=probPer, col=meanA, group=meanA) + 
  geom_line() +
  geom_point() +
  labs(x=expression(paste("log"[10],"(d)")), 
       y=expression(paste("P(N"[i],">0 | N"[i,0],">0)")),
       col="a")

ggsave(paste0("../figures/case2.pdf"), probExc, 
       width=3, height = 3, device = "pdf")  

## Summarize simulated data and add predictions and plot -----
Sims %>% #selection of Sims
  filter(dispType == "exponentialD", meanA<1) %>%
  select(all_of(c("rep", "n", "meanA", "d", "propPatchesN", "vary", "cvA", "k"))) %>%
  group_by(n, meanA, d, vary, cvA, k) %>%
  summarise(meanProb = mean(propPatchesN, na.rm = T), 
            sdProb = sd(propPatchesN, na.rm = T)) %>%
  ggplot() + 
  scale_alpha_manual(values=c(1, 0.4)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=log10(d), y=meanProb, col=meanA, 
                 alpha=as.factor(k)), cex=2) +
  geom_errorbar(aes(x=log10(d), ymin=meanProb-sdProb, 
                    ymax=meanProb+sdProb, 
                    col=meanA, alpha=as.factor(k)), width=0.1) +
  geom_line(data=Predictions, 
            aes(x=log10(d), y=prob, col=meanA, lty=as_factor(meanA)), 
            show.legend = F) +
  facet_grid(cols=vars(cvA), rows=vars(vary), 
             labeller = label_bquote(cols=paste("cv(", a[ij],")=", .(cvA)),
                          rows=paste("vary=", .(vary)))) +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy", col=expression(paste(bar(a))), 
       alpha="k")

ggsave(paste0("../figures/feasExponential.pdf"), width=5, height = 3, 
       device = "pdf")  

## Probability for extinction w/o disp. (theory and sims) --------
Sims %>%
  select(n, d, p, meanA, rep, NHat, vary, k, cvA) %>%
  filter(d==min(d), meanA>0) %>%
  select(-d) %>%
  mutate(fractionExct = map_dbl(NHat, ~.x %>%
                          filter(sp==1, density<extinctionThreshold)%>%
                          nrow)/p) %>%
  group_by(n, p, meanA, cvA, k, vary) %>%
  summarise(meanProb = mean(fractionExct),
            sdProb = sd(fractionExct)) %>%
  left_join(Predictions %>% filter(d==min(d))%>%select(-d), 
            by=c("meanA", "p", "n")) %>%
  ggplot() + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_alpha_manual(values=c(1, 0.4)) + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=meanA, y=meanProb, alpha=as.factor(k))) + 
  geom_errorbar(aes(x=meanA, ymin=meanProb-sdProb, width=0.05,
                    ymax=meanProb+sdProb, alpha=as.factor(k))) +
  aes(x=meanA, y=probN0iExt, group=as.factor(k), 
      alpha=as.factor(k)) + #,,  
  geom_line(show.legend = F) + 
  #scale_color_gradient(low = "yellow", high = "red") +
  facet_grid(cols=vars(cvA), rows=vars(vary), labeller = label_both) +
  labs(x="a", 
       y="Probability that sp. i \n gets excluded w/o disp.", 
       alpha="k")

ggsave(paste0("../figures/PNi0Excl.pdf"), width=6, height = 3, 
       device = "pdf")  

# Leftovers -----
## showcase accuracy of mean r -----
ggplot(IntPredSimple %>% mutate(meanA2 = meanA) %>%
         unite("it", rep, meanA2)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=m, y=meanRPer, col=meanA, pch=as.factor(k)), alpha=0.5) + 
  aes(x=m, y=meanRPerPredicted, col=meanA, lty=as_factor(it)) + #,,  
  geom_line(show.legend = F) + 
  facet_grid(cols=vars(cvA), rows=vars(vary), labeller = label_both) +
  labs(x="nr of persisting species, m", 
       y="mean r of persisting species", 
       col="a", pch="k (equivalence)")

ggsave(paste0("../figures/rm.pdf"), width=6, height = 3, 
       device = "pdf")
#good prediction, except when m<=2 and a is strong. 
#Changing the uniform approximation to a triangular one
#does not seem to improve things (but double-check). So I suspect there is something 
#fundamentally off at m is small. 

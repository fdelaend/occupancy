
# SIMULATIONS ----
## Simulations ------
data <- expand_grid(n = c(4, 6, 8), meanA = c(0.2, 0.4, 0.6, 0.8), #, 6; 0.4, 
                    d = seq(-6,-4, length.out=6), # 6
                    sdA = 0, p = 100, rep = c(1:3)) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  #Make parameters
  mutate(d = 10^d) %>%
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% #make symmetric and ditch diagonal; set self to 1
      make_symmetric() %>% set_diagonal(d=1))) %>% 
  mutate(A = pmap(., make_block_diagonal)) %>% #make A spatial
  mutate(R = pmap(., make_R_spatial)) %>% #make spatial R (note that the emigration is subtracted later so this is the real local R)
  mutate(D = pmap(., make_D)) %>% #make dispersal matrix
  mutate(N0 = map(R, ~ .x*0+extinctionThreshold)) %>% #set initial conditions
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

## Plot simulations  ----
ggplot(data) + 
  scale_color_gradient(low = "yellow", high = "red") +
  scale_linetype_discrete(rep("solid", 100)) +
  theme_bw() +
  aes(x=log10(d), y=propPatchesN, col=meanA) + 
  geom_point() +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy, simulated", col="a") +
  #geom_line(aes(x=log10(d), y=feasibility, col=meanA,
  #              group=interaction(meanA, rep))) + 
  facet_grid(cols=vars(n))

ggsave(paste0("../figures/feas.pdf"), width=6, height = 3, 
       device = "pdf")

## Recover the distribution of the growth rates and plot --------
Rs     <- data %>% filter(d==min(d), meanA==0.8, rep==1) %>% select(n, R, NHat) %>% 
  mutate(NHat = map2(R, NHat, ~ .y %>% mutate(R=.x))) %>%
  #filter(n==6) %>%
  select(-R) %>%
  unnest(NHat) %>%
  select(-density, -sp)

ggplot(Rs) + 
  theme_bw() +
  aes(R, col=as.factor(location)) + 
  geom_density(show.legend=F) + 
  geom_density(aes(R), show.legend=F, lwd=2, col="black") + 
  facet_grid(cols=vars(n))

pdfRs  <- density(Rs$R, from=0)
meanR  <- sum(pdfRs$x*pdfRs$y)/sum(pdfRs$y)#grant mean of R

# INTERMEDIATE PREDICTIONS ----
## Select data w/o dispersal and make predictions ------
dataNoDisp <- data %>%
  filter(d==min(d), meanA>0) %>%
  select(n, meanA, p, rep, summaryM, NTotalK) %>%
  mutate(NTotalK = map_dbl(NTotalK, ~mean(.x$NTotalK))) %>% #mean across sp
  mutate(summaryM = pmap(., function(summaryM, meanA, n, ...) { #add predictions
    summaryM %>%
    mutate(m=as.numeric(as.character(m))) %>%
    rowwise() %>%
    mutate(fractionPatchesPredicted = get_fraction_m(meanA=meanA, m=m, n=n)) %>% #predicted total density in a patch of m persisting sp
    ungroup() %>%
    mutate(meanRPerPredicted = get_RMeanM(a=meanA, m=m, n=n, r=meanR),#predicted mean r of persisting sp
           meanRExcPredicted = (-meanRPerPredicted*m+n*meanR)/(n-m),#predicted mean r of excluded sp
           NTotalPredicted = get_N_total(meanA=meanA, n=m, r=meanRPerPredicted))})) %>% #predicted total density for a patch with m species
  mutate(NTotalKPredicted = p/n*map_dbl(summaryM, ~sum(.x$fractionPatchesPredicted * .x$NTotalPredicted)))

## showcase accuracy of NtotalK ----
ggplot(dataNoDisp) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=NTotalK, y=NTotalKPredicted, col=meanA) + 
  geom_point() + 
  geom_abline(slope=1, intercept=0) + 
  labs(x=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"simulated")), 
       y=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"analytical")), 
       col="a") 

ggsave(paste0("../figures/Nk.pdf"), width=3, height = 2, 
       device = "pdf")
  
## showcase accuracy of f(m) -----
dataNoDispSimple <- dataNoDisp %>%
  select(all_of(c("meanA", "n", "p", "rep", "summaryM"))) %>%
  unnest(summaryM) %>%
  mutate(fractionPatches = nrPatches/p) 

ggplot(dataNoDispSimple %>% mutate(meanA2 = meanA) %>%
         unite("it", rep, meanA2)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  geom_point(aes(x=m, y=fractionPatches, col=meanA), alpha=0.5) + 
  aes(x=m, y=fractionPatchesPredicted, col=meanA, lty=as_factor(it)) + #,,  
  geom_line(show.legend = F) + 
  facet_grid(cols=vars(n)) +
  labs(x="nr of persisting species, m", 
       y="fraction of patches occupied, f(m)", 
       col="a")

ggsave(paste0("../figures/fm.pdf"), width=6, height = 3, 
       device = "pdf")

## showcase accuracy of mean r -----
ggplot(dataNoDispSimple %>% mutate(meanA2 = meanA) %>%
         unite("it", rep, meanA2)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  geom_point(aes(x=m, y=meanRPer, col=meanA), alpha=0.5) + 
  aes(x=m, y=meanRPerPredicted, col=meanA, lty=as_factor(it)) + #,,  
  geom_line(show.legend = F) + 
  facet_grid(cols=vars(n)) +
  labs(x="nr of persisting species, m", 
       y="mean r of persisting species", 
       col="a")

ggsave(paste0("../figures/rm.pdf"), width=6, height = 3, 
       device = "pdf")
#good prediction, except when m<=2 and a is strong. 
#Changing the uniform approximation to a triangular one
#does not seem to improve things (but double-check). So I suspect there is something 
#fundamentally off at m is small. 

# ANALYTICAL PREDICTIONS -----
## Predictions ----
prob <- dataNoDisp %>%
  filter(rep==1) %>%
  mutate(sampleSize = 100) %>%
  mutate(samples = pmap(., function(summaryM, sampleSize, meanA, NTotalKPredicted, p, ...){
    tibble(m = sample(x=summaryM$m, size=sampleSize, prob = summaryM$fractionPatchesPredicted, replace=T)) %>% 
      left_join(summaryM, by="m", multiple="all") %>% #get variables that match the sampled m
      select(all_of(c("m", "NTotalPredicted", "meanRPerPredicted", "meanRExcPredicted"))) %>%
      rowwise() %>%       #sample from distribution of R such that IGR<0 or >0 
      mutate(riExc=sample_ri(samplesize=1, PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="above"),
             riPer=sample_ri(samplesize=1, PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="below")) %>%
      ungroup() %>%
      mutate(N0i=get_N0i(a=meanA, n=m, r=meanRPerPredicted, ri=riPer), 
             N1iExc=get_N1iExc(NTotalK=NTotalKPredicted, ri=riExc, a=meanA, NTotalPredicted),
             Eps = NTotalKPredicted/p*(p-1)-(p-1)*N0i,
             EpsN0i = Eps/N0i)})) %>%
  mutate(means = map(samples, ~ .x %>% group_by(m) %>% #mean across samples with the same m to mimick mean across species 
                       summarize(meanEpsN0i = mean(EpsN0i, na.rm=T),
                                 meanN1iExc=mean(N1iExc, na.rm=T),
                                 meanEps = mean(Eps, na.rm=T)))) %>%
  mutate(samples = map2(samples, means, ~.x %>% #add means to samples
                          left_join(.y, by="m"))) %>%
  expand_grid(d = seq(-6,-4, length.out=6)) %>%
  mutate(d=10^d) %>%
  mutate(samples = pmap(., function(samples, d, meanA, n, ...){ samples %>%
      mutate(N1iPer = get_N1iPer(a=meanA, n=n, m=m, EpsN0i, #N1i when i persists w/o disp.
                                 meanN1iExc, meanEpsN0i), 
             NiExc = d*N1iExc, #Ni when i is excluded w/o disp.
             NiPer = N0i+d*N1iPer)})) %>% #Ni when i persists w/o disp.
  mutate(probExc = map_dbl(samples, ~sum(.x$NiExc>extinctionThreshold)/length(.x$NiExc)), #Compute probabilities that > threshold
         probPer = map_dbl(samples, ~sum(.x$NiPer>extinctionThreshold)/length(.x$NiPer))) %>%
  mutate(summaryM = map2(summaryM, n, ~.x %>%
                             mutate(probN0iExt = fractionPatchesPredicted*(1-m/.y)))) %>%#proba that Ni is absent from patches with m sp (w/o disp)
  mutate(probN0iExt = map_dbl(summaryM, ~sum(.x$probN0iExt)), #overall proba across all patches
         probN0iPer = 1-probN0iExt,
         prob = (probExc*probN0iExt + probPer*probN0iPer)^n) #grant prob

## Plot predictions of main prob ----
ggplot(prob) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=log10(d), y=prob, col=meanA) + 
  geom_point() + 
  facet_grid(cols=vars(n)) +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy, analytical", col="a")

ggsave(paste0("../figures/Analytical.pdf"), width=6, height = 3, 
       device = "pdf")  

## Plot predictions of prob, conditional on i persisting w/o disp. ----
ggplot(prob) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=log10(d), y=probPer, col=meanA) + 
  geom_point() + 
  facet_grid(cols=vars(n)) +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy, analytical", col="a")

ggsave(paste0("../figures/probPer.pdf"), width=6, height = 3, 
       device = "pdf")  

## Add simulated data to the predictions and plot -----
dataSel <- data %>% #selection of data
  select(all_of(c("n", "meanA", "d", "propPatchesN"))) %>%
  left_join(prob, by=c("n", "meanA", "d"))

ggplot(dataSel %>% mutate(meanA2 = meanA) %>%
         unite("it", rep, meanA2)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  geom_point(aes(x=log10(d), y=propPatchesN, col=meanA)) + 
  geom_line(aes(x=log10(d), y=prob, col=meanA, lty=as_factor(it)), 
            show.legend = F) +
  facet_grid(cols=vars(n)) +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy", col="a")

ggsave(paste0("../figures/feasPred.pdf"), width=6, height = 3, 
       device = "pdf")  

## Probability for extinction w/o disp. (theory and sims) --------
dataNoDisp <- data %>%
  select(n, d, p, meanA, rep, NHat) %>%
  filter(d==min(d), meanA>0) %>%
  select(-d) %>%
  mutate(fractionExct = map_dbl(NHat, ~.x %>%
                          filter(sp==1, density<extinctionThreshold)%>%
                          nrow)/p) %>%
  left_join(prob%>%filter(d==min(d))%>%select(-d), 
            by=c("meanA", "rep", "p", "n"))

ggplot(dataNoDisp) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  geom_point(aes(x=meanA, y=fractionExct), alpha=0.5) + 
  aes(x=meanA, y=probN0iExt, lty=as.factor(rep)) + #,,  
  geom_line(show.legend = F) + 
  #scale_color_gradient(low = "yellow", high = "red") +
  facet_grid(cols=vars(n)) +
  labs(x="a", 
       y="Probability that sp. i \n gets excluded w/o disp.", col="a")

ggsave(paste0("../figures/PNi0Excl.pdf"), width=6, height = 3, 
       device = "pdf")  


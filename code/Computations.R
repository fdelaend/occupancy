
# CALCULATIONS ------
data <- expand_grid(n = c(4, 6, 8), meanA = c(0.2, 0.4, 0.6, 0.8), #, 6; 0.4, 
                    d = seq(-8,-4, length.out=2), # 6
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

## fit distribution to r
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
  #geom_line(aes(x=log10(d), y=feasibility, col=meanA,
  #              group=interaction(meanA, rep))) + 
  facet_grid(cols=vars(n))

ggsave(paste0("../figures/feas.pdf"), width=6, height = 2, 
       device = "pdf")

## get data w/o dispersal and make predictions ------
dataNoDisp <- data %>%
  filter(d==min(d), meanA>0) %>%
  select(n, meanA, p, rep, summaryM, NTotalK) %>%
  mutate(NTotalK = map_dbl(NTotalK, ~mean(.x$NTotalK))) %>% #mean across sp
  mutate(summaryM = pmap(., function(summaryM, meanA, n, ...) {
    summaryM %>%
    mutate(m=as.numeric(as.character(m))) %>%
    rowwise() %>%
    mutate(fractionPatchesPredicted = get_fraction_m(meanA=meanA, m=m, n=n)) %>% #predicted total density in a patch of m persisting sp
    ungroup() %>%
    mutate(meanRPerPredicted = get_RMeanM(a=meanA, m=m, n=n, r=1.02),#predicted mean r of persisting sp
           meanRExcPredicted = (-meanRPerPredicted*m+n*1.02)/(n-m),#predicted mean r of excluded sp
           NTotalPredicted = get_N_total(meanA=meanA, n=m, r=meanRPerPredicted))})) %>% #total density for a patch with m species
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
       col="a") #+
  #facet_grid(cols=vars(n))

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

## Analytical predictions -----
Prob <- dataNoDisp %>%
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
  #NOW group by m to compute the means of N1iExc and EpsN0i
  expand_grid(d=seq(-6,-4, length.out=6)) %>%
  mutate(d=10^d) %>%
  mutate(samples = pmap(., function(d, p, NTotalKPredicted, meanA, samples, ...){
    samples %>% mutate(NiExc=d*N1iExc,#Ni for when i is extinct w/o dispersal
                       #meanN1iExc = mean(N1iExc), NOT CORRECT, is across sp..
                       #meanEps = mean(Eps), #NOT CORRECT, is across sp..
                       #meanEpsN0i = mean(EpsN0i), NOT CORRECT, is across sp..
                       N1iPer = 5)})) %>% 
  mutate(ProbNi = map_dbl(samples, ~sum(.x$NiExc>extinctionThreshold)/length(.x$NiExc))) %>%#Compute probability that Ni>0
  mutate(summaryM = map2(summaryM, n, ~.x %>%
                             mutate(probNi0Ext = fractionPatchesPredicted*(1-m/.y)))) %>%#proba that Ni is absent from patches with m sp (w/o disp)
  mutate(ProbNi0Ext = map_dbl(summaryM, ~sum(.x$probNi0Ext))) #overall proba across all patches
  
ggplot(Prob1) + 
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red") +
  aes(x=log10(d), y=ProbNi, col=meanA) + 
  geom_point() + 
  facet_grid(cols=vars(n)) +
  labs(x=expression(paste("log"[10],"(d)")), 
       y="Patch occupancy sp i, \n exclusion w/o dispersal", col="a")

ggsave(paste0("../figures/PNiExcl.pdf"), width=6, height = 3, 
       device = "pdf")  

#Now probability for extinction w/o disp. (theory and sims)
dataNoDisp <- data %>%
  select(n, d, p, meanA, rep, NHat) %>%
  filter(d==min(d), meanA>0) %>%
  select(-d) %>%
  mutate(fractionExct = map_dbl(NHat, ~.x %>%
                          filter(sp==1, density<extinctionThreshold)%>%
                          nrow)/p) %>%
  left_join(Prob1%>%filter(d==min(d))%>%select(-d), 
            by=c("meanA", "rep", "p", "n"))
  
ggplot(dataNoDisp) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_gradient(low = "yellow", high = "red") +
  geom_point(aes(x=meanA, y=fractionExct, ), alpha=0.5) + 
  aes(x=meanA, y=ProbNi0Ext, lty=as.factor(rep)) + #,,  
  geom_line(show.legend = F) + 
  #scale_color_gradient(low = "yellow", high = "red") +
  facet_grid(cols=vars(n)) +
  labs(x="a", 
       y="Fraction of patches where i \n gets excluded w/o disp.", col="a")

ggsave(paste0("../figures/PNi0Excl.pdf"), width=6, height = 3, 
       device = "pdf")  

#  mutate(r = map(ri, ~mean(.x)),  #compute mean and mean of inverse
#  rInv = map(ri, ~mean(1/.x))) %>%
#  mutate(m = map2(mProbs, sampleSize, ~sample(x=c(1:length(.x)), size=.y, prob=.x, replace=T))) %>% #sample m's
#  mutate(NiWeak = pmap(., get_density_weak)) %>% #compute Ni for weak interactions
#mutate(Ni = pmap(., function(meanA, NiWeak, NiStrong,...){NiWeak*(meanA<0.5)+NiStrong*(meanA>=0.5)})) %>%#get proper Ni according to meanA
#  mutate(ProbNi0Zero = pmap_dbl(., function(n, mProbs, ...){
#sum(mProbs[1:(n-1)]*(1-c(1:(n-1))/n)) })) 



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

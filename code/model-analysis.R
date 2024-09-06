
source("tools-simulations.R")
source("tools-other.R")
# Read simulation output and put them into a single object
Sims <- tibble(filenr = c(1:100)) |>
  mutate(data = map(filenr, ~read_simulations(paste("../simulated-data/",.x,"data.rds", sep="")))) |>
  unnest(data) |>
  select(!filenr) |>
  mutate(rep = as.numeric(rep))
# Save the object
#saveRDS(Sims, file=paste("../simulated-data/all-data.RDS",sep=""))
Sims <- readRDS(file="../simulated-data/all-data.RDS")
## Recover the distribution of the growth rates --------
Rs     <- Sims |> 
  #Only keep a factorial combo, since the same Rs are being sampled everywhere
  filter(meanA==0.2, d==min(d), vary==0, k==1, cvA==0, p==100, 
         rep==1) |> 
  select(R) |>
  unnest(R)
pdfRs  <- density(Rs$R, from=0)
meanR  <- sum(pdfRs$x*pdfRs$y)/sum(pdfRs$y)#grant mean of R

# INTERMEDIATE PREDICTIONS: ----
## Make the predictions of fm and NTotalK -----
IntPred <- Sims |> #Start with simulations to obtain the predictor values (to make sure we're predicting for exactly the same parameter values) 
  # Select Simulations w/o dispersal (select regularD although the same for both dispersal types),  
  # stable matrices, diffuse competition, and 
  # regional equivalence, and make predictions
  filter(meanA<1, d==min(d), vary==0, k==1, cvA==0, 
         dispType=="regularD") |>
  # *Measure* total density across all patches, mean across species (NTotalK)
  mutate(NTotalK = map_dbl(NTotalK, ~mean(.x$NTotalK))) |> #mean regional density across sp
  # *Predict*, for each level of local richness m: 
  # -fm: expected fraction of patches with richness=m (fmPredicted),
  # -mean r of persisting sp (meanRPerPredicted)
  # -mean r of excluded sp (meanRExcPredicted) 
  # -total density across all patches, mean across species (NTotalPredicted)
  (\(x) mutate(x, summaryM = pmap(x, \(summaryM, meanA, n, ...)
                                  summaryM |>
                                    mutate(m=as.numeric(as.character(m))) |>
                                    rowwise() |>
                                    mutate(fmPredicted = get_fraction_m(meanA=meanA, m=m, n=n)) |> #predicted total density in a patch of m persisting sp
                                    ungroup() |>
                                    mutate(meanRPerPredicted = get_RMeanM(a=meanA, m=m, n=n, r=meanR),#predicted mean r of persisting sp
                                           meanRExcPredicted = (-meanRPerPredicted*m+n*meanR)/(n-m),#predicted mean r of excluded sp
                                           NTotalPredicted = get_N_total(meanA=meanA, n=m, r=meanRPerPredicted)))))() |> #predicted total density for a patch with m species
  # Predict total regional density (same for all species)
  mutate(NTotalKPredicted = p/n*map_dbl(summaryM, ~ (.x |> 
                                          mutate(NperM = fmPredicted * NTotalPredicted)|>
                                          summarise(sum(NperM)))[[1]])) #Total regional density for a species

## Unnest and only keep those variables that are for now relevant ---
SelIntPred <- IntPred |>
  unnest(summaryM) |>
  select(all_of(c("n", "meanA", "d", "vary", "k", "cvA", "dispType",
                  "p", "rep", "m", "fmPredicted")))

## Add predictions of fm to simulations -----
IntPredsAndSims <- Sims |>
  unnest(summaryM) |>
  mutate(m=as.numeric(as.character(m))) |>
  mutate(fractionPatches = nrPatches/p) |> #get simulated fm
  select(all_of(c("d", "dispType", "meanA", "n", "p", "rep", "m", "fractionPatches", 
                  "cvA", "vary", "k"))) |> #only keep relevant variables
  left_join(SelIntPred, by=c("d", "meanA", "n", "p", "rep", "m", "cvA", "vary", "k", "dispType"),
            multiple = "all") |> #join with predictions; 
  #makes sure that predictions are NA where dispersal is > min(d),
  #or where there is no regional equivalence
  #compute summary stats of predicted and simulated fm 
  #(mean and sd)
  summarise(meanProb = mean(fractionPatches),
            sdProb = sd(fractionPatches),
            meanProbPred = mean(fmPredicted),
            sdProbPred = sd(fmPredicted), 
            .by=c("d", "dispType", "meanA", "n", "p", "m", "cvA", "vary", "k"))

## plot simulations vs. predictions of fm; only take k==1 and vary==0 -----
IntPredsAndSims |> 
  mutate(d=round(log10(d),1)) |>
  filter(k==1, dispType=="regularD", 
         d %in% c(-6, -5.2, -4.4), vary==0, p==100) |> 
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
       y="fraction of patches with m species, fm", 
       col=expression(paste(bar(a))))

ggsave(paste0("../figures/fm.pdf"), width=4.5, height = 3, 
       device = "pdf")

## Plot predicted vs. simulated NtotalK ----
NtotalK <- ggplot(IntPred |>
                    summarise(meanNTotalK = mean(NTotalK),
                              sdNTotalK = sd(NTotalK),
                              meanNTotalKPredicted = mean(NTotalKPredicted),
                              sdNTotalKPredicted = sd(NTotalKPredicted),
                              .by = c(meanA, p))) + 
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=meanA, y=log10(meanNTotalK), group=as.factor(p)) + 
  geom_point() + 
  geom_errorbar(aes(x=meanA, ymin=log10(meanNTotalK-sdNTotalK), 
                    ymax=log10(meanNTotalK+sdNTotalK), 
                    width = 0.05)) + 
  geom_line(aes(x=meanA, y=log10(meanNTotalKPredicted), 
                lty=as.factor(p)), show.legend = F) +
  labs(x="interaction strength, a", 
       y=expression(paste("log"[10]," of total density, ", Sigma[{k}],"N"[{"0,i"}]^{(k)}))) 
#labs(x=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"simulated")), 
#     y=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"analytical"))) +
#facet_grid(cols=vars(cvA), rows=vars(vary), labeller = label_both)

ggsave(paste0("../figures/Nk.pdf"), width=4, height = 3, 
       device = "pdf")

# PREDICTIONS OF PATCH OCCUPANCY -----
## Make the predictions for when all assumptions are met
Predictions <- IntPred |>
  filter(rep==1) |> #Select 1st replicate
  select(-c("rep", "k", "cvA", "vary", "d")) |> #Remove irrelevant variables
  mutate(sampleSize = 1000) |> #set sample size for probability calculations
  (\(x) mutate(x, samples = pmap(x, sample_random)))() |>
  mutate(means = map(samples, ~ .x |>  
                       summarize(meanN1iExc=mean(N1iExc, na.rm=T), #mean across patches of Ni1 in case of exclusion w/o dispersal. Because we sample 1 species per patch, this is the same as taking a mean across species
                                 rho = mean(rhoi, na.rm=T),
                                 .by = m))) |> #same type of mean, but now of rhoi 
  mutate(samples = map2(samples, means, ~.x |> #add means to samples
                          left_join(.y, by = join_by(m)))) |>
  expand_grid(d = seq(-6,-4, length.out=60)) |> #Create a gradient of dispersal values
  mutate(d=10^d) |>
  (\(x) mutate(x, samples = pmap(x, sample_random_Ni)))() |>
  mutate(probExc = map_dbl(samples, ~sum(.x$NiExc>extinctionThreshold)/length(.x$NiExc)), #Compute probabilities that > threshold
         probPer = map_dbl(samples, ~sum(.x$NiPer>extinctionThreshold)/length(.x$NiPer))) |>
  mutate(summaryM = map2(summaryM, n, ~.x |>
                             mutate(probN0iExt = fmPredicted*(1-m/.y)))) |> #proba that Ni is absent from patches with m sp (w/o disp)
  mutate(probN0iExt = map_dbl(summaryM, ~sum(.x$probN0iExt)), #overall proba across all patches
         probN0iPer = 1-probN0iExt,
         prob = (probExc*probN0iExt + probPer*probN0iPer)^n) |> #grant prob
  mutate(fm = map(summaryM, ~ .x |> select(all_of(c("m", "nrPatches"))) |>
                    mutate(fm = nrPatches/sum(nrPatches)) |>
                    select(all_of(c("m", "fm"))))) |>
  (\(x) mutate(x, A = pmap(x, make_A)))() |>
  mutate(Xi = map_dbl(A, ~feasibility(.x))) |>
  (\(x) mutate(x, prob2 = pmap_dbl(x, get_patch_occupancy)))() |>
  select(-summaryM)

## Illustrate predictions of NtotalK and negative IGR ----
NegIGR <- Predictions |> 
  filter(d==min(d)) |>
  select(all_of(c("meanA", "samples", "p"))) |>
  unnest(samples) |>
  ggplot() +
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=NegIGR, col=meanA, lty=as.factor(p), 
      group=interaction(meanA, as.factor(p))) + 
  geom_density(show.legend = F) + 
  labs(x="negative invasion growth rate", 
       y="probability density", col="a")

## Plot predictions of prob, conditional on i persisting/excluded w/o disp. ----
probExc <- ggplot(Predictions) + 
  theme_bw() +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=log10(d), y=probExc, col=meanA, lty=as.factor(p), 
      group=interaction(meanA, as.factor(p))) + 
  geom_line() +
  #geom_point() +
  labs(x=expression(paste("log"[10],"(D)")), 
       y=expression(paste("P(N"[i],">0 | N"[i,0],"<0)")),
       col="a", lty="nr of patches")

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
       y=expression(paste("P(N"[i],">0 | N"[i0],">0)")),
       col="a")

ggsave(paste0("../figures/case2.pdf"), probExc, 
       width=3, height = 3, device = "pdf")  

## Summarize simulated data and add predictions and plot -----
SimsSum <- Sims |> 
  select(all_of(c("p", "rep", "n", "meanA", "d", 
                  "propPatchesN", "vary", "cvA", "k", "dispType"))) %>%
  group_by(n, meanA, d, vary, cvA, k, p, dispType) %>%
  summarise(meanProb = mean(propPatchesN, na.rm = T), 
            sdProb = sd(propPatchesN, na.rm = T))

#Plot for when all assumptions are met
ggplot(SimsSum |>
         filter(dispType == "regularD", meanA<1, k==1, cvA==0, vary==0)) + 
  scale_alpha_manual(values=c(0.3, 0.6, 1)) + 
  theme_bw() +
  scale_linetype_manual(values=rep("solid", 1000)) + 
  scale_color_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=log10(d), y=meanProb, col=meanA), cex=2) +
  geom_errorbar(aes(x=log10(d), ymin=meanProb-sdProb, 
                    ymax=meanProb+sdProb, 
                    col=meanA), width=0.1) +
  geom_line(data=Predictions, 
            aes(x=log10(d), y=prob2, group=interaction(meanA, as.factor(p)), col=meanA), 
            show.legend = F) +
  facet_grid(.~p, 
             labeller = label_bquote(cols=paste(.(p)," patches"))) +
  labs(x=expression(paste("log"[10],"(D)")), 
       y="Patch occupancy", col="a", 
       alpha="k")

ggsave(paste0("../figures/feasIdeal.pdf"), width=6, height = 2, 
       device = "pdf")  

#Plot for when not met
ggplot(SimsSum |> filter(dispType == "exponentialD", vary == 0)) + 
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
            aes(x=log10(d), y=prob2, col=meanA, group=as_factor(meanA)), 
            show.legend = F) +
  facet_grid(cvA ~ p, labeller = label_bquote(rows=paste("cv(", a[ij],")=", .(cvA)),
                                     cols=paste(.(p)," patches"))) +
  labs(x=expression(paste("log"[10],"(D)")), 
       y="Patch occupancy", col=expression(paste(bar(a))), 
       alpha="k")

ggsave(paste0("../figures/feasExponential.pdf"), width=6, height = 3, 
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

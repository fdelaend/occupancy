
source("tools-simulations.R")
source("tools-other.R")
# Read simulation output and put them into a single object
Sims <- tibble(filenr = c(1:100)) |>
  mutate(data = map(filenr, ~read_rds(paste("../simulated-data/",.x,"data.rds", sep="")))) |>
  unnest(data) |>
  select(!filenr) |>
  mutate(rep = as.numeric(rep))
# Save the object
#saveRDS(Sims, file=paste("../simulated-data/all-data.RDS",sep=""))
Sims <- readRDS(file="../simulated-data/all-data.RDS")
## Recover the distribution of the growth rates --------
Rs     <- Sims |> 
  #Only keep a factorial combo, since the same Rs are being sampled everywhere
  filter(meanA==0.2, d==min(d), vary==0, k==1, cvA==0, p==80, 
         rep==1) |> 
  select(R) |>
  unnest(R)
pdfRs  <- density(Rs$R, from=0)
meanR  <- sum(pdfRs$x*pdfRs$y)/sum(pdfRs$y)#grant mean of R

# PLOTS TO GET INTUITION -----
Sims |> 
  filter(meanA %in% c(0.2, 0.8), vary == 0, log10(d) < -4,
         k == 1, cvA == 0.2, p == 40, 
         dispType == "regularD") |>
  select(!R & !NHat & !NTotalK) |>
  unnest(summaryM) |>
  mutate(m = as.numeric(as.character(m))) |>
  mutate(fraction = nrPatches / p) |>
  summarise(median_fraction = median(fraction), 
            min_fraction = quantile(fraction, 0.25),
            max_fraction = quantile(fraction, 0.75),
            .by = c(m, d, meanA)) |>
  ggplot() +
  scale_color_viridis_d(option="plasma", end=0.9) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x=log10(d), y=median_fraction, col=as.factor(m)) +
  geom_line() + 
  #geom_errorbar(aes(x=log10(d), ymin=min_fraction, 
  #                  ymax=max_fraction, col=as.factor(m), 
  #                  width = 0.1)) +
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste(bar(a),"=", .(meanA)))) +
  labs(y = "fraction", col = "local richness", x = "log10 of dispersal rate")

# PREDICTIONS of fm and NTotalK, for diffuse competition ----
## Predictions of fm -----
Predictions_fm <- expand_grid(n = 6, m = c(1:6), iteration = c(1:10),
                   meanA = seq(0, 0.8, 0.05)) |>
  mutate(meanA = if_else(meanA==1, meanA+1e-5, meanA)) |> #Avoid matrices with det()=0
  (\(x) mutate(x, fmPredicted = pmap_dbl(x, get_fraction_m)))() |>
  #fm isn't calculated correctly for m = 1, so get it here as 1 - sum of fm for 2 to 6
  mutate(fmPredicted = if_else(m==1, abs(1-sum(fmPredicted, na.rm=T)), fmPredicted),
         .by = c(n, meanA, iteration))  |>
  summarise(fmPredicted = mean(fmPredicted), 
            .by = c(n, m, meanA))

#saveRDS(Predictions_fm, file=paste("../simulated-data/Predictions_fm.RDS",sep=""))
Predictions_fm <- readRDS(file="../simulated-data/Predictions_fm.RDS")
  
## Predictions of NTotalK ----
Predictions_NK <- Predictions_fm |> 
  # Predict, for each level of local richness m: 
  # -mean r of persisting sp (meanRPerPredicted)
  # -mean r of excluded sp (meanRExcPredicted) 
  # -total density across all patches, mean across species (NTotalPredicted)
  mutate(meanRPerPredicted = get_RMeanM(a=meanA, m=m, n=n, r=meanR),#predicted mean r of persisting sp
         meanRExcPredicted = (-meanRPerPredicted*m+n*meanR)/(n-m),#predicted mean r of excluded sp
         NTotalPredicted = get_N_total(meanA=meanA, n=m, r=meanRPerPredicted)) |>
  # Predict total regional density (same for all species)
  expand_grid(p = seq(10, 100, 1)) |>
  mutate(NperM = p/n*fmPredicted * NTotalPredicted) |>
  mutate(NTotalKPredicted = sum(NperM), 
         .by = c(n, meanA, p)) 

NtotalK <- ggplot(Predictions_NK |>
                    filter(meanA<1) |>
                    summarise(meanNTotalKPredicted = mean(NTotalKPredicted),
                              .by = c(meanA, p))) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  geom_tile(aes(x=meanA, fill=log10(meanNTotalKPredicted), 
                y=p)) +
  labs(x="interaction strength, a", 
       fill=expression(paste("log"[10],"(", Sigma[{k}],"N"[{"0,i"}]^{(k)}, ")"))) +
  theme(legend.position="bottom")
#labs(x=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"simulated")), 
#     y=expression(paste(Sigma[{k}],"N"[{"0,i"}]^{(k)},~"analytical"))) +
#facet_grid(cols=vars(cvA), rows=vars(vary), labeller = label_both)

ggsave(paste0("../figures/Nk.pdf"), width=3, height = 3, 
       device = "pdf")

# PREDICTIONS OF PATCH OCCUPANCY -----
## Make the predictions for when all assumptions are met
Predictions <- Predictions_NK |>
  filter(meanA > 0) |> #Otherwise error when computing sample_ri
  nest_by(n, meanA, p, NTotalKPredicted) |>
  ungroup() |>
  mutate(sampleSize = 1000) |> #set sample size for probability calculations
  (\(x) mutate(x, samples = pmap(x, sample_random)))() |>
  mutate(samples = map(samples, ~ .x |>  
                       mutate(meanN1Exc=mean(N1iExc, na.rm=T), #mean across patches of Ni1 in case of exclusion w/o dispersal. Because we sample 1 species per patch, this is the same as taking a mean across species
                              meanrho = mean(rhoi, na.rm=T), #same type of mean, but now of rhoi
                              .by = m))) |>  
  expand_grid(d = seq(-6,-0, length.out=10)) |> #Create a gradient of dispersal values
  mutate(d=10^d) |>
  (\(x) mutate(x, samples = pmap(x, sample_random_Ni)))() |> #Sample random values for Ni (density of i with dispersal)
  mutate(probExc = map_dbl(samples, ~sum(.x$NiExc>extinctionThreshold)/length(.x$NiExc)), #Compute probabilities that > threshold
         probPer = map_dbl(samples, ~sum(.x$NiPer>extinctionThreshold)/length(.x$NiPer))) |>
  (\(x) mutate(x, A = pmap(x, make_A)))() |>
  mutate(Xi = map_dbl(A, ~feasibility(.x))) |>
  (\(x) mutate(x, prob2 = pmap_dbl(x, get_patch_occupancy)))() |>
  select(!data & !A)

saveRDS(Predictions, file=paste("../simulated-data/Predictions.RDS",sep=""))

## Illustrate predictions of NtotalK and negative IGR ----
NegIGR <- Predictions |> 
  filter(d==min(d)) |>
  select(all_of(c("meanA", "samples", "p"))) |>
  unnest(samples) |>
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=NegIGR, col=meanA, group = meanA) + 
  geom_density() + 
  labs(x="exclusion rate", 
       y="probability density", col="a") +
  theme(legend.position="bottom")

## Plot predictions of prob, conditional on i excluded w/o disp. ----
probExc <- ggplot(Predictions |> 
                    filter(meanA %in% c(0.1, 0.8), log10(d) < -4)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  aes(x=log10(d), y=p, fill=probExc) + 
  geom_tile() +
  #geom_point() +
  labs(x=expression(paste("log"[10],"(D)")), 
       fill=expression(paste("P(N"[i],">0 | N"[i0],"=0)")),
       y="nr of patches") + 
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste("a = ",.(meanA)))) +
  theme(legend.position="bottom")

case1 <- grid.arrange(NtotalK, NegIGR, probExc, 
                      layout_matrix = rbind(c(1, 2), c(3, 3)), 
                      widths=c(1,1)) 

ggsave(paste0("../figures/case1.pdf"), case1, 
       width=5, height = 7.5, device = "pdf")  

## Plot Case 2
probExc <- ggplot(Predictions |> 
                    filter(meanA %in% c(0.1, 0.8), log10(d) < -4)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  aes(x=log10(d), y=p, fill=probPer) + 
  geom_tile() +
  labs(x=expression(paste("log"[10],"(D)")), 
       y="nr of patches", 
       fill=expression(paste("P(N"[i],">0 | N"[0i],">0)")),
       col="a") +
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste("a = ",.(meanA)))) +
  theme(legend.position="bottom")

ggsave(paste0("../figures/case2.pdf"), probExc, 
       width=3, height = 3, device = "pdf")  

## Summarize simulated data to add to predictions -----
SimsSum <- Sims |> 
  select(all_of(c("n", "meanA", "d", "vary", "k", 
                  "cvA", "p", "dispType", "propPatchesN"))) |>
  summarise(meanProb = mean(propPatchesN, na.rm = T), 
            sdProb = sd(propPatchesN, na.rm = T), 
            .by = c(n, meanA, d, vary, k, cvA, p, dispType)) |>
  mutate(k = if_else(k==1, "Regional equivalence", "Regional dominance"))

#Plot for when all assumptions are met, but allowing d to be large
ggplot(SimsSum |>
         filter(dispType == "regularD", 
                meanA<1, k=="Regional equivalence", cvA==0, vary==0)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  geom_tile(data=Predictions |> filter(meanA %in% c(0.2, 0.4, 0.8)), 
              aes(x=log10(d), y=p, fill = prob2), show.legend = F) +
  geom_point(aes(x=log10(d), y=p, fill=meanProb), col = "white", 
             pch = 21, cex=3) +
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste("a = ", .(meanA)))) +
  labs(x=expression(paste("log"[10],"(D)")), 
       y="nr. of patches", fill="Patch occupancy")

ggsave(paste0("../figures/feasIdeal.pdf"), width=6, height = 2, 
       device = "pdf")  

#Plot for when not met
ggplot(SimsSum |> filter(dispType == "exponentialD", 
                         vary > 0, cvA == 0.2)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=log10(d), y=p, fill = meanProb), pch = 21, cex=5) +
  facet_grid(k~meanA, 
             labeller = label_bquote(cols = paste(bar(a), " = ", .(meanA)))) +
  labs(x = expression(paste("log"[10],"(D)")), 
       y = "nr of patches, p", fill= "Patch occupancy")

ggsave(paste0("../figures/feasExponential.pdf"), width=6, height = 3, 
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

source("tools-simulations.R")
source("tools-other.R")
# Compile all simulation output from cluster into a single object ---------
# or load the compiled data if it already exists
output_file <- "../simulated-data/simulated-data/all-data.RDS"

if (file.exists(output_file)) {
  Sims <- readRDS(file=output_file)
} else {
  Sims <- tibble(filenr = c(1:100)) |>
    mutate(data = map(filenr, ~read_rds(paste("../simulated-data/simulated-data/",.x,"data.rds", sep="")))) |>
    unnest(data) |>
    select(!filenr) |>
    mutate(rep = as.numeric(rep))

  # Save the object
  saveRDS(Sims, file=output_file)
  }

# Recover the distribution of the growth rates --------
Rs     <- Sims |> 
  #Only keep a single combination, since the same Rs are being sampled everywhere
  filter(meanA==0.5, d==min(d), vary==0, k==1, cvA==1e-2, p==50, 
         rep==1) |> 
  select(R) |>
  unnest(R)
pdfRs  <- density(Rs$R, from=0)
meanR  <- sum(pdfRs$x*pdfRs$y)/sum(pdfRs$y)#grant mean of R


# Predict f(m) (eq.14 in SI) for diffuse competition ----
# or load if it already exists
output_file <- "../simulated-data/simulated-data/Predictions_fm.RDS"

if (file.exists(output_file)) {
  Predictions_fm <- readRDS(file=output_file)
} else {
  Predictions_fm <- expand_grid(n = 10, m = c(1:10), 
                                iteration = c(1:10),
                                meanA = seq(0.1, 0.8, 0.05)) |>
    mutate(meanA = if_else(meanA==1, meanA+1e-5, meanA)) |> #Avoid matrices with det()=0
    (\(x) mutate(x, fmPredicted = pmap_dbl(x, get_fraction_m)))() |>
    #fm isn't calculated correctly for m = 1, so get it here as 1 - sum of fm for 2 to 6
    mutate(fmPredicted = if_else(m==1, abs(1-sum(fmPredicted, na.rm=T)), fmPredicted),
           .by = c(n, meanA, iteration))  |>
    summarise(fmPredicted = mean(fmPredicted), 
              .by = c(n, m, meanA))
  # Save the object
  saveRDS(Predictions_fm, file=output_file)
}

# Predict regional total of density of a species (eq.13 in SI) for diffuse competition ----
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
                y=p), colour = "NA") +
  labs(x="interaction strength, a", y = "nr of patches, p", 
       fill=expression(paste("log"[10],"(", Sigma[{k}],hat(N)[{"0,i"}]^{(k)}, ")"))) +
  theme(legend.position="bottom")

# Predict probability to persist when excluded w/o dispersal -----
# or load if it already exists
output_file <- "../simulated-data/simulated-data/Predictions.RDS"

if (file.exists(output_file)) {
  Predictions <- readRDS(file=output_file)
} else {
  Predictions <- Predictions_NK |>
    filter(meanA > 0) |> #Otherwise error when computing sample_ri
    nest_by(n, meanA, p, NTotalKPredicted) |>
    ungroup() |>
    mutate(sampleSize = 100) |> #set sample size for probability calculations. 
    (\(x) mutate(x, samples = pmap(x, sample_random)))() |>
    mutate(samples = map(samples, ~ .x |>  
                           mutate(meanN1Exc=mean(N1iExc, na.rm=T), #mean across patches of Ni1 in case of exclusion w/o dispersal. Because we sample 1 species per patch, this is the same as taking a mean across species
                                  meanrho = mean(rhoi, na.rm=T), #same type of mean, but now of rhoi
                                  .by = m))) |>  
    expand_grid(d = seq(-6,-2, length.out=20)) |> #Create a gradient of dispersal values. 
    mutate(d=10^d) |>
    (\(x) mutate(x, samples = pmap(x, sample_random_Ni)))() |> #Sample random values for Ni (density of i with dispersal)
    mutate(probExc = map_dbl(samples, ~sum(.x$NiExc>extinctionThreshold)/length(.x$NiExc)), #Compute probabilities that > threshold
           probPer = map_dbl(samples, ~sum(.x$NiPer>extinctionThreshold)/length(.x$NiPer))) |>
    (\(x) mutate(x, A = pmap(x, make_A)))() |>
    mutate(Xi = map_dbl(A, ~feasibility(.x))) |>
    (\(x) mutate(x, prob2 = pmap_dbl(x, get_patch_occupancy)))() |>
    select(!data & !A)
  saveRDS(Predictions, file=output_file)
}

# Plot this probability ----
probExc <- ggplot(Predictions |> 
                    filter(meanA %in% c(0.2, 0.5), log10(d) < -4)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  aes(x=log10(d), y=p, fill=probExc) + 
  geom_tile(color = NA) +
  #geom_point() +
  labs(x=expression(paste("dispersal rate, log"[10],"(d)")), 
       fill=expression(paste("P(",hat(N)[i],">",10^{-3}," | ",hat(N)["0i"],"=0)")),
       y="nr of patches, p") + 
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste("a = ",.(meanA)))) +
  theme(legend.position="bottom")

# Get out the predictions of invasion growth rate (IGR)
Predictions_IGR <- Predictions |>
  filter(d==min(d), p == 100) |>
  select(all_of(c("meanA", "samples", "p"))) |>
  unnest(samples)

# Illustrate predictions of negative invasion growth rate ----
NegIGR <- Predictions_IGR |> 
  ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_viridis_c(option="plasma", end=0.9) +
  aes(x=NegIGR, col=meanA, group = meanA) + 
  geom_density() + 
  labs(x="exclusion rate", 
       y="probability density", col="a") +
  theme(legend.position="bottom")

# Fig 2 ----
Fig2 <- (NtotalK | NegIGR) / probExc + plot_annotation(tag_levels = "A") 

ggsave(paste0("../figures/Fig2.pdf"), Fig2, 
       width=4.2, height = 6.3, device = "pdf")  

# Evidence that persistence with disperal is almost guaranteed ----
# when persistence happenss without dispersal.
probPer <- ggplot(Predictions |> 
                    filter(meanA %in% c(0.2, 0.5), log10(d) < -4)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9, 
                       n.breaks = 2) +
  aes(x=log10(d), y=p, fill=probPer) + 
  geom_tile() +
  labs(x=expression(paste("dispersal rate, log"[10],"(d)")), 
       y="nr of patches", 
       fill=expression(paste("P(N"[i],">0 | N"["0i"]," > 0)")),
       col="a") +
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste("a = ",.(meanA)))) +
  theme(legend.position="bottom")

ggsave(paste0("../figures/FigS3.pdf"), probPer, 
       width=4, height = 3, device = "pdf")  

# Summarize simulated data to add to predictions -----
SimsSum <- Sims |> 
  select(all_of(c("n", "meanA", "d", "vary", "k", 
                  "cvA", "p", "dispType", "propPatchesN"))) |>
  summarise(meanProb = mean(propPatchesN, na.rm = T), 
            sdProb = sd(propPatchesN, na.rm = T), 
            .by = c(n, meanA, d, vary, k, cvA, p, dispType)) |>
  mutate(k = if_else(k==1, "Equivalence", "Inequivalence"))

#Fig3A: Plot for when all assumptions are met -----
Fig3A <- ggplot(SimsSum |>
         filter(dispType == "regularD", meanA < 1, 
                d < 1e-2, 
                k=="Equivalence", cvA==0.01, vary==0)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  geom_tile(data=Predictions |> filter(meanA %in% c(0.2, 0.5), 
                                       p <= 51), 
              aes(x=log10(d), y=p, fill = prob2), show.legend = F) +
  geom_point(aes(x=log10(d), y=p, fill=meanProb), col = "white", 
             shape = 22, cex=2) +
  facet_grid(.~meanA, 
             labeller = label_bquote(cols=paste("a = ", .(meanA)))) +
  labs(x=expression(paste("dispersal rate, log"[10],"(d)")), 
       y="nr. of patches, p", fill="patch occupancy")


#Fig 3B: Plot for when assumptions are not met ------
Fig3B <- ggplot(SimsSum |> filter(dispType == "exponentialD", 
                         vary > 0, cvA == 0.2)) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(colour = "snow2")) +
  scale_fill_viridis_c(option="plasma", end=0.9) +
  geom_point(aes(x=log10(d), y=p, fill=meanProb), col = "white", 
             shape = 22, cex=3) +
#  geom_tile(aes(x=log10(d), y=p, fill = meanProb)) +#, pch = 21, cex=3
  facet_grid(k~meanA, 
             labeller = label_bquote(cols = paste(bar(a), " = ", .(meanA)))) +
  labs(x = expression(paste("dispersal rate, log"[10],"(d)")), 
       y = "nr of patches, p", fill= "patch occupancy")

#Fig3 ---------
Fig3 <- Fig3A / Fig3B +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  plot_annotation(tag_levels = "A")

ggsave(paste0("../figures/Fig3.pdf"), 
       width=5, height = 6.5, device = "pdf")  

# Showcase accuracy of mean r -----
Mean_predictions <- Predictions_NK |>
  filter(p==50, meanA %in% c(0.2, 0.5)) |>
  summarise(pred = mean(meanRPerPredicted, na.rm=T),
            .by = c(m, meanA))

Sims |> 
  filter(dispType == "regularD", meanA < 1, d==10^-6, 
         k==1, cvA==0.01, vary==0, p==50) |>
  select(meanA, summaryM) |>
  unnest(summaryM) |>
  mutate(m = as.numeric(as.character(m))) |>
  ggplot() + 
  theme_bw() +
  scale_colour_viridis_d(option="plasma", end=0.9) +
  aes(x = m, y = meanRPer, col = as.factor(meanA)) + 
  geom_point() + 
  geom_line(data = Mean_predictions, 
             aes(x = m, y = pred, col = as.factor(meanA))) + 
  labs(x = "nr of persisting species, m", 
       y = "mean growth rate", 
       col = "a")

ggsave(paste0("../figures/FigS2.pdf"), width=4, height = 3, 
       device = "pdf")

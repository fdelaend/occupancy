
source("tools-simulations.R")
source("tools-other.R")

# Combine all simulation output from cluster into a single object ----
Sims <- tibble(filenr = c(1:100)) |>
  mutate(data = map(filenr, ~read_rds(paste("../simulated-data/simulated-data-daphnids/",.x,"data.rds", sep="")))) |>
  unnest(data) |>
  select(!filenr) |>
  mutate(rep = as.numeric(rep))
# Save the object
saveRDS(Sims, file=paste("../simulated-data/simulated-data-daphnids/all-data.RDS",sep=""))
# Or read it if already done ----
Sims <- readRDS(file="../simulated-data/simulated-data-daphnids/all-data.RDS")

# Summarize simulated data -----
SimsSum <- Sims |> 
  select(all_of(c("n", "meanA", "d", "vary", "k", 
                  "cvA", "p", "dispType", "summaryM", "propPatchesN", "propPatches2"))) |>
  summarise(meanProb = mean(propPatchesN, na.rm = T), 
            sdProb = sd(propPatchesN, na.rm = T), 
            meanProb2 = mean(propPatches2, na.rm = T), 
            sdProb2 = sd(propPatches2, na.rm = T), 
            .by = c(n, meanA, d, vary, k, cvA, p, dispType))# |>
  #mutate(k = if_else(k==1, "Equivalence", "Inequivalence"))

#Plot -----
# choice of d range based on Dubart. Max of about 10e-4. 
# Only when d gets at the higher end: p has negative effect.
ggplot(SimsSum |>
         filter(meanA == 1.2, cvA == 0.01, vary == 0.1)) + 
  theme_bw() +
  scale_colour_gradient(low = "grey85", high = "black") +
  theme(panel.grid = element_blank()) +
  geom_line(aes(x=p, y=meanProb2, col = log10(d),
                group = log10(d))) +
  facet_grid(k~dispType) +
  labs(col=expression(paste("dispersal rate, log"[10],"(d)")), 
       x="nr. of patches, p", y="patch occupancy")

ggsave(paste0("../figures/patch-occupancy-daphnids.pdf"), 
       width=6, height = 3, device = "pdf")  

ggplot(SimsSum |>
         filter(vary > 0, cvA == 0.2, 
                log10(d) > -6, log10(d)< -4)) + 
  theme_bw() +
  scale_colour_gradient(low = "grey85", high = "black") +
  theme(panel.grid = element_blank()) +
  geom_line(aes(x=log10(d), y=meanProb2, col = p,
                group = p)) +
  facet_grid(k~dispType) +
  labs(x=expression(paste("dispersal rate, log"[10],"(d)")), 
       col="nr. of patches, p", y="patch occupancy")


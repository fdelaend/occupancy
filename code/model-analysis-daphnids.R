
source("tools-simulations.R")
source("tools-other.R")
# Compile all simulation output from cluster into a single object ---------
# or load the compiled data if it already exists
output_file <- "../simulated-data/simulated-data-daphnids/all-data.RDS"

if (file.exists(output_file)) {
  Sims <- readRDS(file=output_file)
} else {
  Sims <- tibble(filenr = c(1:100)) |>
    mutate(data = map(filenr, ~read_rds(paste("../simulated-data/simulated-data-daphnids/",.x,"data.rds", sep="")))) |>
    unnest(data) |>
    select(!filenr) |>
    mutate(rep = as.numeric(rep))
  saveRDS(Sims, file=output_file)
}

# Summarize simulated data -----
SimsSum <- Sims |> 
  select(all_of(c("n", "meanA", "d", "vary", "k", 
                  "cvA", "p", "dispType", "summaryM", 
                  "propPatchesN", "propPatches2"))) |>
  summarise(meanProb = mean(propPatchesN, na.rm = T), 
            sdProb = sd(propPatchesN, na.rm = T), 
            meanProb2 = mean(propPatches2, na.rm = T), 
            sdProb2 = sd(propPatches2, na.rm = T), 
            .by = c(n, meanA, d, vary, k, cvA, p, dispType)) |>
  mutate(k = if_else(k==1, "Equivalence", paste("Inequivalence, k = ", k)))

# Only select realistic parameter settings ; ------
# These are those for which meanProb < 0.01 across all values of dispersal 
SimsSum_realistic <- SimsSum |>
  group_by(d, meanA, vary, k, cvA, dispType) |>
  filter(all(meanProb < 0.01)) |>
  ungroup()

#Fig 7: Plot vs. p -----
ggplot(SimsSum_realistic |> 
         filter(dispType == "exponentialD")) + # "exponentialD"
  theme_bw() +
  scale_colour_viridis_c(option="plasma", end=0.9) +
  theme(panel.grid = element_blank()) +
  aes(x=p, y=meanProb2, col = meanA) +
  geom_point() +
  facet_grid(k~cvA, 
             labeller = label_bquote(cols = paste(cv(a), " = ", .(cvA)))) + #k, cvA, dispType
  labs(col=expression(paste(bar(a))), 
       x="nr. of patches, p", y="patch occupancy of pairs when n=3")

ggsave(paste0("../figures/Fig7.pdf"), 
       width=6, height = 6, device = "pdf")  

#Fig 8: Plot vs. d -----
ggplot(SimsSum_realistic |> 
         filter(dispType == "exponentialD")) + # "exponentialD"
  theme_bw() +
  scale_colour_viridis_c(option="plasma", end=0.9) +
  theme(panel.grid = element_blank()) +
  aes(x=as.factor(log10(d)), y=meanProb2, col = meanA,
      group = interaction(log10(d), meanA)) +
  geom_boxplot() +
  facet_grid(k~cvA, 
             labeller = label_bquote(cols = paste(cv(a), " = ", .(cvA)))) + #k, cvA, dispType
  labs(col=expression(paste(bar(a))), 
       x=expression(paste("dispersal rate, log"[10],"(d)")), 
       y="patch occupancy of pairs when n=3")

ggsave(paste0("../figures/Fig8.pdf"), 
       width=6, height = 6, device = "pdf")  

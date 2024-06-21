library(tidyverse)


calculateCoexistenceProb <- function(dat) {
  sumDat <- dat |>
    summarise(count = sum(count), .by = c(sample, taxon)) |>
    filter(count > 0) |>
    summarise(n_taxa = sum(count > 0), .by = sample)
  tibble(num_persisting_taxa = sort(unique(sumDat$n_taxa))) |>
    mutate(prob_coex = map_dbl(num_persisting_taxa, \(m) {
      nrow(filter(sumDat, n_taxa == m)) / nrow(sumDat)
    }))
}


dat <- read_rds("../data/all-data.rds") |>
  # Unselect overly fine-grained taxonomic levels:
  select(!genus & !family & !order) |>
  # Unite remaining taxonomic level names into single taxon string:
  unite(col = "taxon", !file & !sample & !count)

coexTable <- dat |>
  nest(data = !file) |>
  mutate(probCoex = map(data, calculateCoexistenceProb)) |>
  select(!data) |>
  unnest(probCoex)

coexTable |>
  ggplot(aes(x = num_persisting_taxa, y = prob_coex)) +
  geom_point(colour = "steelblue", alpha = 0.5) +
  facet_wrap(~ file, scales = "free") +
  labs(x = "Number of coexisting taxa", y = "Probability of coexistence") +
  theme_bw() +
  theme(panel.grid = element_blank())

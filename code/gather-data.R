library(tidyverse)
library(fs)


# Each data file has a similar structure, so they can be cleaned in a uniform way:
cleanData <- function(dataFromFile) {
  dataFromFile |>
    # Sometimes the sample ID is purely numeric; convert to character string:
    mutate(Sample = as.character(Sample)) |>
    # Rename column to lowercase:
    rename(sample = Sample) |>
    # Tidy up the data (it will have 3 columns: sample, taxon, and count):
    pivot_longer(!sample, names_to = "taxon", values_to = "count") |>
    # Make sure "count" is an integer:
    mutate(count = as.integer(count))
}


tibble(file = Sys.glob("../data/*.tsv")) |>
  # Load all data files:
  mutate(data = map(file, compose(as_tibble, read_tsv))) |>
  mutate(file = path_file(file)) |>
  # Clean up each data file:
  mutate(data = map(data, cleanData)) |>
  unnest(data) |>
  # Some samples have taxon = "Unclassified"; remove these:
  filter(taxon != "Unclassified") |>
  # Clean up taxon names by removing "p__" for phylum, "o__" for order, etc.:
  mutate(taxon = str_remove_all(taxon, ".__")) |>
  # Separate out taxonomic classification in separate columns:
  separate_wider_delim(cols = taxon, delim = ";",
                       names = c("domain","phylum","class","order","family","genus")) |>
  write_rds("../data/all-data.rds", compress = "xz")

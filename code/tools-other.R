library(tidyverse)
library(brms)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(pracma)
library(gridExtra)
library(feasoverlap)
library(patchwork)
library(reticulate) #to run python code from R
use_python("/usr/local/bin/python3.10") #load your preferred python version
source_python("MVN.py")
library(grateful)
library(xtable)
library(geodist)
library(sf)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggspatial)
library(cowplot)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#

#' Equilibrium density of species i without dispersal.
#'
#' @param a Interspecific competition coefficient.
#' @param n Number of species in the local community.
#' @param r Mean intrinsic growth rate of non-focal species.
#' @param ri Intrinsic growth rate of the focal species.
#'
#' @return Equilibrium density of species i (Eq. 21).
get_N0i <- function(a=0.5, n=10, r=1, ri=0.1){
  (a*(n - 1)*r - ri*(a*(n - 2) + 1))/((a - 1)*(a*(n - 1) + 1))
}

#' Total equilibrium abundance without dispersal.
#'
#' @param meanA Mean interspecific interaction strength.
#' @param d Self-regulation term.
#' @param n Number of species.
#' @param r Mean intrinsic growth rate.
#'
#' @return Total equilibrium abundance (Eq. 22).
get_N_total <- function(meanA=0.2, d=1, n=10, r=1){ 
  n*r/(d+meanA*(n-1))
}

#' Density contribution of unit dispersal for excluded species.
#'
#' @param NTotalK Total abundance across other patches.
#' @param ri Intrinsic growth rate of focal species.
#' @param a Interaction strength.
#' @param SumN0j Total abundance of resident species.
#'
#' @return Density contribution per unit dispersal (Eq. 7).
get_N1iExc <- function(NTotalK=10, ri=1, a=0.5, SumN0j){
  NTotalK/(a*SumN0j - ri)
}

#' Density contribution of unit dispersal for persisting species.
#'
#' @param a Interaction strength.
#' @param n Number of species in the regional pool.
#' @param m Number of species in the focal patch.
#' @param meanrho Mean rho across species.
#' @param rhoi Rho of focal species.
#' @param p Number of patches.
#' @param meanN1Exc Mean N1 for excluded species.
#'
#' @return Density contribution per unit dispersal (Eq. 29).
get_N1iPer <- function(a=0.5, n=4, m=1, meanrho, rhoi, p, meanN1Exc){
  (a*((1-a)*meanN1Exc*(n-m)+(m-1)*meanrho+(2-m)*rhoi) + (p-1)*(1-a) - rhoi)/
    ((1-a)*(a*(m-1)+1))
}

#' Fraction of patches with m species (no dispersal).
#'
#' @param meanA Interaction strength.
#' @param m Number of species in a patch.
#' @param n Regional species pool size.
#'
#' @return Fraction of patches with m species (Eq. 14).
get_fraction_m <- function(meanA=0.5, m=2, n=3, ...){
  A   <- diag(n) + meanA
  A   <- set_diagonal(A, 1)
  Am <- A %*% diag(c(rep(1,m),rep(0,n-m)))
  diag(Am) <- c(rep(1,m), rep(-1,n-m)) 
  B  <- diag(c(rep(1, n)))
  feas1combo <- 2^(n)*calculate_omega_constraint(A=Am, B=B)^(n)
  ifelse(m==1, NA, choose(n, m) * feas1combo)
}

#' Distribution of species richness across patches.
#'
#' @param n Number of species.
#' @param meanA Interaction strength.
#'
#' @return Vector of f(m) for m = 1,...,n.
make_distribution <- function(n, meanA=0.5){
  probs <- NULL
  for (m in c(1:n))
    probs <- c(probs, get_fraction_m(meanA=meanA, m=m, n=n))
  probs
}

## Truncate a discretized probability density.
#'
#' @param pdfFitted Output of density().
#' @param q Quantile threshold.
#' @param ditch Direction of truncation ("up" or "down").
#'
#' @return Truncated density.
trunc_dist <- function(pdfFitted, q = 0.5, ditch = "up") {
  cumprob <- cumsum(pdfFitted$y) / sum(pdfFitted$y)
  Qs <- if (ditch == "down") which(cumprob >= q) else which(cumprob <= q)
  list(x=pdfFitted$x[Qs], y=pdfFitted$y[Qs])
}

## Mean of a truncated density.
#'
#' @param pdfFitted Output of trunc_dist().
#' @param q Quantile threshold.
#' @param ditch Direction of truncation.
#'
#' @return Mean of truncated distribution.
get_mean_trunc <- function(pdfFitted, q=0.5, ditch="up"){
  pdfTrunc <- trunc_dist(pdfFitted, q=q, ditch=ditch)
  sum(pdfTrunc$x*pdfTrunc$y)/sum(pdfTrunc$y)
}

#' Mean intrinsic growth rate given local richness.
#'
#' @param a Interaction strength.
#' @param m Local species richness.
#' @param n Regional species pool size.
#' @param r Mean intrinsic growth rate.
#' @param x Scaling factor.
#'
#' @return Mean intrinsic growth rate (Eq. 25).
get_RMeanM <- function(a=0.8, m=1, n=4, r=1.01, x=2){
  ((1 + a*(-1 + m))*n*r*x)/(m*(a*(n + m*(-1 + x) - x) + x)) 
}

## Sample intrinsic growth rates from truncated distribution.
#'
#' @param samplesize Number of samples.
#' @param PDF Density with components x and y.
#' @param cutoff Truncation threshold.
#' @param ditch Direction ("below" or "above").
#'
#' @return Sampled intrinsic growth rate(s).
sample_ri <- function(samplesize=1, PDF, cutoff, ditch="below", ...){
  ditchBelow <- sample(PDF$x[PDF$x>cutoff], size=1, prob=PDF$y[PDF$x>cutoff], replace=T)
  ditchAbove <- sample(PDF$x[PDF$x<cutoff], size=1, prob=PDF$y[PDF$x<cutoff], replace=T)
  (ditch=="below")*ditchBelow + (ditch=="above")*ditchAbove
}

## Patch occupancy prediction.
#'
#' @param data Tibble with f(m).
#' @param probExc Persistence probability given exclusion without dispersal.
#' @param n Regional species pool size.
#' @param Xi Feasibility for n species.
#'
#' @return Patch occupancy probability.
get_patch_occupancy <- function(data, probExc, n, Xi, ...){
  data |>
    filter(m<n) |>
    mutate(term = fmPredicted*probExc^(n-m)) |>
    summarize(sum(term)) |>
    (\(x) x*(1-Xi) + Xi)() |>
    as_vector()
}

## Sample random variables for patch occupancy prediction.
#'
#' @return Tibble of sampled quantities used in analytical predictions.
sample_random <- function(data, sampleSize, meanA, NTotalKPredicted, p, ...){
  tibble(m = sample(x=data$m, size=sampleSize, prob=data$fmPredicted, replace=T)) %>%
    left_join(data, by="m", multiple="all") %>%
    select(all_of(c("m","NTotalPredicted","meanRPerPredicted","meanRExcPredicted"))) %>%
    rowwise() %>%
    mutate(riExc=sample_ri(PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="above"),
           riPer=sample_ri(PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="below")) %>%
    ungroup() %>%
    mutate(N0i=get_N0i(a=meanA, n=m, r=meanRPerPredicted, ri=riPer),
           N1iExc=get_N1iExc(NTotalK=NTotalKPredicted, ri=riExc, a=meanA, NTotalPredicted),
           NegIGR=1/get_N1iExc(NTotalK=1, ri=riExc, a=meanA, NTotalPredicted),
           rhoi = NTotalKPredicted/p*(p-1)/N0i)
}

## Sample equilibrium density with dispersal.
#'
#' @return Tibble with Ni for excluded and persisting cases.
sample_random_Ni <- function(samples, d, meanA, n, p, ...){ 
  samples |>
    mutate(N1iPer = get_N1iPer(a=meanA, n=n, m=m, meanrho=meanrho,
                               rhoi=rhoi, meanN1Exc=meanN1Exc, p=p),
           NiExc = d*N1iExc,
           NiPer = N0i+d*N1iPer)
}
 
## Fit Bayesian GLMM for priority effects.
#'
#' @return Fitted brms model.
fit_priority <- function(data, focal = "magna",
                         adapt_delta = 0.9,
                         warmup = 5000, iter = 7000) {
  form <- paste0(focal, "_summer_alone ~ as.factor(", focal, "_spring)",
                 " + (1|island/pool) + (1|year)")
  brm(as.formula(form),
      data = data,
      family = brms::bernoulli(),
      control = list(adapt_delta = adapt_delta),
      warmup = warmup,
      iter = iter)
}

## Model printing for paper 
#' The ... contain the variables to be selected apart from 'model'
print_model <- function(data, ...){
  test <- data |>
    ungroup() |>
    mutate(results = map(model, ~ tidy(.x))) |>
    select(results, ...) |>
    unnest(results) |>
    as.data.frame() |>
    xtable() |>
    print(type = "latex",
          scientific = TRUE,
          include.rownames = FALSE)
}



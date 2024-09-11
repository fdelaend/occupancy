library(mvtnorm)
library(pracma)
library(gridExtra)
library(feasoverlap)
library(reticulate) #to run python code from R
use_python("/usr/local/bin/python3.10") #load your preferred python version
source_python("MVN.py")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#

get_N_total <- function(meanA=0.2, d=1, n=10, r=1){ #a=comp.strength; d=self limit.
  n*r/(d+meanA*(n-1))
}

#get the density of a species i when there is no dispersal
get_N0i <- function(a=0.5, n=10, r=1, ri=0.1){
  (a*(n - 1)*r - ri*(a*(n - 2) + 1))/((a - 1)*(a*(n - 1) + 1))
}

#get N1i when i is excluded w/o dispersal
#NTotalK = total abundance of i across all patches without the focal patch
#ri = r of i, given that it gets excluded w/o dispersal
#a = interaction strength
#SumN0j = abundance of all species 
get_N1iExc <- function(NTotalK=10, ri=1, a=0.5, SumN0j){
  NTotalK/(a*SumN0j - ri)
}

#get N1i when i persists w/o dispersal
#a = interaction strength
#n and m: nr of species in the regional pool and in the focal patch
#meanN1iExc = mean of N1iExc across species
#rho = mean of rhoi across species
#rhoi = rho for species i, given that it persists w/o dispersal
get_N1iPer <- function(a=0.5, n=4, m=1, rho, rhoi, p, meanN1iExc){
  (a*((1-a)*meanN1iExc*(n-m)+(m-1)*rho+(2-m)*rhoi) + (p-1)*(1-a) - rhoi)/
    ((1-a)*(a*(m-1)+1))
}

#compute fraction of patches with m species in a metacommunity of n species
#without dispersal. For m>1
get_fraction_m <- function(meanA=0.5, m=2, n=3, ...){
  #n x n matrix
  A   <- diag(n) + meanA
  A   <- set_diagonal(A, 1)
  Am <- A %*% diag(c(rep(1,m),rep(0,n-m)))
  #Am2<- diag(c(rep(1,m),rep(0,n-m))) %*% Am
  diag(Am) <- c(rep(1,m), rep(-1,n-m)) 
  #diag(Am2)<- 1
  B        <- diag(c(rep(1, n))) #constraints
  feas1combo <- 2^(n)*calculate_omega_constraint(A=Am, B=B)^(n)
  ifelse(m==1, NA, choose(n, m) * feas1combo)
}

#make a distribution of nr of species in a patch in
#case there is no dispersal
make_distribution <- function(n, meanA=0.5){
  probs <- NULL
  for (m in c(1:n))
  {
    probs <- c(probs, get_fraction_m(meanA=meanA, m=m, n=n))
  }
  probs
}
#truncate a distribution (i.e. cut off fraction below (if ditch="down") 
#or above (ditch="up") a certain quantile q)
#pdfFitted = output of density()
#q = quantile, where 1 = 100th quantile, 0.5 = 50th etc..
trunc_dist <- function(pdfFitted, q=0.5, ditch="up"){
  cumprob <- cumsum(pdfFitted$y)/sum(pdfFitted$y) #discretize to cum. proba.
  #identify x to be ditched
  if(ditch=="down") {Qs<-which(cumprob>=q)} else {Qs<-which(cumprob<=q)}
  list(x=pdfFitted$x[Qs], y=pdfFitted$y[Qs])
}

#get mean of a truncated pdf
#pdfFitted = output of trunc_dist()
#q = quantile, where 1 = 100th quantile, 0.5 = 50th etc..
get_mean_trunc <- function(pdfFitted, q=0.5, ditch="up"){
  pdfTrunc <- trunc_dist(pdfFitted, q=q, ditch=ditch)
  sum(pdfTrunc$x*pdfTrunc$y)/sum(pdfTrunc$y)
}

#mean R of the persisting sp
#x = factor to divide a*Nt by in order to get the mean
get_RMeanM <- function(a=0.8, m=1, n=4, r=1.01, x=2.2){
  ((1 + a*(-1 + m))*n*r*x)/(m*(a*(n + m*(-1 + x) - x) + x)) 
} 

#sample ri's from a truncated dist
# PDF = a pdf, with $x the values and $y the density
# cutoff = threshold for truncation
# ditch = whether to ditch everything below or above
sample_ri <- function(samplesize=1, PDF, cutoff, ditch="below", ...){
  ditchBelow <- sample(x = PDF$x[which(PDF$x>cutoff)], size=1, 
                       prob = PDF$y[which(PDF$x>cutoff)], replace = T)
  ditchAbove <- sample(x = PDF$x[which(PDF$x<cutoff)], size=1, 
                       prob = PDF$y[which(PDF$x<cutoff)], replace = T)
  (ditch=="below")*ditchBelow + (ditch=="above")*ditchAbove
}

# Get patch occupancy
# Input:
# -fm: fraction of patches with m coexisting species
# -m: nr of coexisting species, in same order as fm
# -probExc: probability for a species to persist if it does not persist without dispersal
# -n: total number (i.e. number of regionally available) species
# -Xi: feasibility for n species
get_patch_occupancy <- function(fm, probExc, n, Xi, ...){
  fm |>
    filter(m<n) |>
    mutate(term = fm*probExc^(n-m)) |>
    summarize(sum(term)) |>
    (\(x) x*(1-Xi) + Xi)() |>
    #(\(x) x + Xi)() |>
    as_vector()
}

#Sample the random variables needed to predict patch occupancy
#Input:
# -summaryM: a tibble with the number of species in a patch (m),
# the total biomass (NTotal), given m, 
# the mean growth rate of the persisting species (meanRPer),
# the predicted fraction of patches in which m species persist
#Output:
# - samples: a tibble with samples for 
# m, riExc (growth rate of sp excluded w/o disp), 
# riPer (growth rate of sp persisting w/o disp),
# N0i (density w/o disp of a sp persisting w/o disp)
# N1iExc (contribution of a unit dispersal to density of a sp excluded w/o disp)
# NegIGR (negative invasion growth rate of sp excluded w/o disp)
# rhoi (ratio of local vs. mean regional abundance)
sample_random <- function(summaryM, sampleSize, meanA, NTotalKPredicted, p, ...){
  tibble(m = sample(x=summaryM$m, size=sampleSize, prob = summaryM$fmPredicted, replace=T)) %>% #sample m according to f(m). Every sample is a hypothetical patch      
    left_join(summaryM, by="m", multiple="all") %>% #get variables that match the sampled m
    select(all_of(c("m", "NTotalPredicted", "meanRPerPredicted", "meanRExcPredicted"))) %>%
    rowwise() %>%       
    #sample 1 growth rate per hypo patch, from distribution of R such that IGR<0 (for riExc) or >0 (for riPer) 
    mutate(riExc=sample_ri(samplesize=1, PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="above"),
           riPer=sample_ri(samplesize=1, PDF=pdfRs, cutoff=meanA*NTotalPredicted, ditch="below")) %>%
    ungroup() %>%
    mutate(N0i=get_N0i(a=meanA, n=m, r=meanRPerPredicted, ri=riPer), #predict density of persisting sp in absence of dispersal
           N1iExc=get_N1iExc(NTotalK=NTotalKPredicted, ri=riExc, a=meanA, NTotalPredicted), #density contribution per unit of dispersal, in case of exclusion w/o dispersal
           NegIGR=1/get_N1iExc(NTotalK=1, ri=riExc, a=meanA, NTotalPredicted), #negative IGR
           rhoi = NTotalKPredicted/p*(p-1)/N0i)}

# Samples random values for Ni (density of i with dispersal)
# Input:
# - samples: output from sample_random
# - d (dispersal rate), meanA (mean comp. interaction strength),
# - n (nr of sp), p (nr of patches)
sample_random_Ni <- function(samples, d, meanA, n, p, ...){ 
  samples |>
    mutate(N1iPer = get_N1iPer(a=meanA, n=n, m=m, rho=rho, #N1i when i persists w/o disp.
                               rhoi=rhoi, meanN1iExc=meanN1iExc, p=p), 
           NiExc = d*N1iExc, #Ni when i is excluded w/o disp.
           NiPer = N0i+d*N1iPer)}
 
# Read the simulation results, 
# and ditch the variables that take up lots of memory and won't be used
# and summarize the density data
# (This last step is a workaround until R gets an update on the cluster;
# This summary is currently not done correctly bc of an older R version on the cluster)
read_simulations <- function(file){
  readRDS(file) |>
    select(!summaryM) |>
    #1/ Summarize the simulated data: per m, compute the nr of patches and total biomass of an average patch
    (\(x) mutate(x, summaryM = pmap(x, \(NHat, R, n,...) 
                                    NHat |>
                                      mutate(present = density>extinctionThreshold,
                                             R=R) |>
                                      summarize(m = as_factor(sum(present)),
                                                NTotal = sum(density),
                                                meanRPer = sum(R*present)/sum(present), 
                                                .by = location) |>
                                      mutate(m = fct_expand(m, as.character(c(1:n)))) |>
                                      group_by(m, .drop=F) |>
                                      summarize(nrPatches = n(),#nr of patches with m sp.
                                                NTotal = mean(NTotal),#total biomass in a patch with m sp.
                                                meanRPer = mean(meanRPer)))))()  |>#mean r of persisting sp.
    #2/proportion of patches in which all n species persist
    mutate(propPatchesN = 1/p*map2_dbl(summaryM, n, ~ (.x |> filter(m==.y))$nrPatches)) |>
    #3/total density across all patches of a species
    mutate(NTotalK = map(NHat, ~ .x |> 
                           summarize(NTotalK = sum(density), .by = sp))) |>
    select(!A & !D & !distances & !N0 & !coords) #ditch all unneeded data
}


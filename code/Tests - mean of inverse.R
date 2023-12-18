#a = off-diagonal element of A
get_diagonal_random <- function(a=0.01, n=10){
  (a*(n-2)+1)/((1-a)*(a*(n-1)+1))
}

get_offdiagonal_random <- function(a=0.1, n=10){
  a/((a-1)*(a*(n-1)+1))
}

data <- expand_grid(n = c(5, 10, 20, 40, 80, 160), meanA = 0, 
                    sdA = c(0.05, 0.1, 0.2, 0.4),
                    rep=c(1:10)) %>% 
  #Make parameters
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% 
      #make symmetric and ditch diagonal; set diagonal to 1
      set_diagonal(d=1))) %>% #make_symmetric() %>% 
  mutate(AInv = map(A, ~solve(.x)), #compute inverse
         meanDiagInv = map_dbl(AInv, ~mean(diag(.x))), #mean of diagonals
         meanOffDiagInv = map2_dbl(AInv, n, ~sum(.x*(1-diag(.y)))/(.y^2-.y)), #mean of off-diagonals
         meanDiagInvPred = get_diagonal_random(a=sdA, n=n),
         meanOffDiagInvPred = get_offdiagonal_random(a=sdA, n=n),
         feas = map2_dbl(AInv, n, ~prod((.x%*%rep(1,.y))>0)),#feasibility for r=1 for all i
         Ntotal = map2_dbl(AInv, n, ~sum((.x%*%rep(1,.y)))),
         NtotalPred = get_N_total(meanA=meanA, d=1, n=n, r=1)) 

ggplot(data %>% filter(feas==1)) + 
  aes(x=log10(n), y=Ntotal, pch=as.factor(sdA)) + 
  geom_point() + 
  geom_point(aes(x=log10(n), y=NtotalPred), col="red") +
  labs(x="log(n)", 
       y="N total")

ggplot(data %>% filter(feas==1)) + 
  aes(x=meanDiagInv, y=meanDiagInvPred, pch=as.factor(sdA)) + 
  geom_point() + 
  geom_point(aes(x=meanOffDiagInv, pch=as.factor(sdA), y=meanOffDiagInvPred), 
             col="red") +
  labs(x="obs", 
       y="pred")



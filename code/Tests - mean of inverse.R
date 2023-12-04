
data <- expand_grid(n = c(5, 10, 20, 40, 80, 160), meanA = 0, sdA = 0.1,
                    rep=c(1:10)) %>% 
  #Make parameters
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% 
      #make symmetric and ditch diagonal; set diagonal to 1
      make_symmetric() %>% set_diagonal(d=1))) %>%
  mutate(AInv = map(A, ~solve(.x)), #compute inverse
         meanDiagAInv = map_dbl(AInv, ~mean(diag(.x))), #mean of diagonals
         meanOffDiagInv = map2_dbl(AInv, n, ~sum(.x*(1-diag(.y)))/(.y^2-.y)), #mean of off-diagonals
         feas = map2_dbl(AInv, n, ~prod((.x%*%rep(1,.y))>0))) #feasibility for r=1 for all i

ggplot(data) + 
  aes(x=log10(n), y=meanDiagAInv) + 
  geom_point() + 
  geom_point(aes(x=log10(n)+0.01, y=meanOffDiagInv), col="red") +
  labs(x="log(n)", 
       y="mean of inverse's diags (black) and off-diags (red)")





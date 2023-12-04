
data <- expand_grid(n = c(5, 10, 20, 40, 80, 160), meanA = c(0.1), sdA = 0.1,
                    rep=c(1:10)) %>% 
  #Make parameters
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% #make symmetric and ditch diagonal; set self to 1
      make_symmetric() %>% set_diagonal(d=1))) %>%
  mutate(AInv = map(A, ~solve(.x)),
         meanDiagAInv = map_dbl(AInv, ~mean(diag(.x))),
         meanOffDiagInv = map2_dbl(AInv, n, ~sum(.x*(1-diag(.y)))/(.y^2-.y)),
         feas = map2_dbl(AInv, n, ~prod((.x%*%rep(1,.y))>0)))

ggplot(data) + 
  aes(x=log10(n), y=meanDiagAInv) + 
  geom_point() + 
  geom_point(aes(x=log10(n)+0.01, y=meanOffDiagInv), col="red")





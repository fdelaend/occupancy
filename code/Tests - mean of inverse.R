
data <- expand_grid(n = c(5, 10, 20, 40, 80), meanA = c(0), sdA = 0.5) %>% 
  #Make parameters
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n) %>% #make symmetric and ditch diagonal; set self to 1
      set_diagonal(d=1))) %>%
  mutate(AInv = map(A, ~solve(.x)),
         meanDiagAInv = map_dbl(AInv, ~mean(diag(.x))),
         meanOffDiagInvTimesN = n*map_dbl(AInv, ~1/(nrow(.x)*(nrow(.x)-1))*sum(.x - diag(x=diag(.x), nrow=nrow(.x)))))

ggplot(data) + 
  aes(x=n, y=meanDiagAInv) + 
  geom_point() + 
  geom_point(aes(x=n, y=meanOffDiagInvTimesN), col="red")





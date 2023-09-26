# Just based on available data
dataNoDisp <- data %>%
  filter(d==min(d), n==3, rep==2) %>%
  select(meanA, NHat) %>%
  unnest(NHat) %>%
  mutate(density=(density>extinctionThreshold)*1) %>%
  select(-X) %>%
  pivot_wider(names_from=sp, values_from=c(R,density)) %>%
  filter(meanA==0.4) %>%
  unite("community", contains("density"))


#Now do properly from scratch
dataNoDisp <- expand_grid(n = 3, meanA = c(0.5),
                    d = 0, sdA = 0, p = 200, rep = c(1)) %>% #nr of species, mean and cv of a, nr of patches in landscape; nr of reps
  mutate(A = pmap(., function(meanA, sdA, n, ...) 
    matrix(data = rnorm(n^2, meanA, sdA), ncol = n))) %>% 
  mutate(A = map(A, ~make_symmetric(.x))) %>% #make symmetric and ditch diagonal
  mutate(A = map(A, ~ set_diagonal(A=.x, d=1))) %>%
  mutate(A = pmap(., make_block_diagonal)) %>% #make A spatial
  mutate(rawR = map2(n, p, ~t(sample_sphere_surface(dim=.x, n = .y, 
                                                    radius = 1, positive = TRUE)))) %>%
  mutate(R = map(rawR, ~c(t(.x)))) %>%
  mutate(D = pmap(., make_D)) %>% #make dispersal matrix
  mutate(N0 = map(R, ~ .x*0+extinctionThreshold)) %>% #set initial conditions
  mutate(NHat = pmap(., get_NHat)) %>%
  mutate(NHat = map2(NHat, R, ~cbind(.x, R=.y))) %>%
  select(meanA, NHat) %>%
  unnest(NHat) %>%
  mutate(density=(density>extinctionThreshold)*1) %>%
  pivot_wider(names_from=sp, values_from=c(R,density)) %>%
  unite("community", contains("density")) %>%
  filter(meanA==0.5)

coms <- unique(dataNoDisp$community)
colr <- cbPalette[c(1:length(coms))]
names(colr) <- coms
s3d <- scatterplot3d::scatterplot3d(x=dataNoDisp$R_1, y=dataNoDisp$R_2, 
                                    z=dataNoDisp$R_3, color=colr[dataNoDisp$community], 
                                    type = "p", pch = 10, cex.symbols=2,
                                    box = F, angle = 1, scale.y=0.11)

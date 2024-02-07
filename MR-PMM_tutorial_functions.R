## SIMULATE FROM MR-PMM


mrpmm_sim <- function(n, m, Sigma_b, Sigma_e, phy) {

  n = n # number of species
  m = m # number of traits
  C <- ape::vcv.phylo(phy, corr = T) # scaled covariance matrix
  I = diag(n) # identity matrix
  u = rep(0,m) # fixed effects (intercepts for each trait)
  
  # simulate phylogenetic random effects
  b = mvrnorm(n = 1, mu = rep(0,n*m), Sigma = kronecker(Sigma_b,C))
  # simulate errors
  e = mvrnorm(1,rep(0,n*m),kronecker(Sigma_e,I))
  # define reps
  rep = length(b)/m
  # generate data frame
  animal <- phy$tip.label
  d <- data.frame(animal)
  
  b.list <- list()
  e.list <- list()
  for (i in 1:m){ # group effects by response trait
    b.list[[i]] <- b[((1+rep*i)-rep):(rep*i)]
    e.list[[i]] <- e[((1+rep*i)-rep):(rep*i)]
  }
  names(b.list) <- c(paste0("b",1:m))
  names(e.list) <- c(paste0("e",1:m))
  
  for (i in 1:m){ # construct response traits from each linear predictor
    d[paste0("y",i)] = u[i] + unlist(b.list[i]) + unlist(e.list[i])
  }
    
  return (d)

}

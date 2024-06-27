## MR-PMM TUTORIAL FUNCTIONS


mrpmm_sim <- function(n, m, Sigma_phy, Sigma_res, phy) {

  n = n # number of species
  m = m # number of traits
  C <- ape::vcv.phylo(phy, corr = T) # scaled covariance matrix
  I = diag(n) # identity matrix
  u = rep(0,m) # fixed effects (intercepts for each trait)
  
  # simulate phylogenetic random effects
  u = mvrnorm(n = 1, mu = rep(0,n*m), Sigma = kronecker(Sigma_phy,C))
  # simulate errors
  e = mvrnorm(1,rep(0,n*m),kronecker(Sigma_res,I))
  # define reps
  rep = length(u)/m
  # generate data frame
  animal <- phy$tip.label
  d <- data.frame(animal)
  
  u.list <- list()
  e.list <- list()
  for (i in 1:m){ # group effects by response trait
    u.list[[i]] <- u[((1+rep*i)-rep):(rep*i)]
    e.list[[i]] <- e[((1+rep*i)-rep):(rep*i)]
  }
  names(u.list) <- c(paste0("u",1:m))
  names(e.list) <- c(paste0("e",1:m))
  
  for (i in 1:m){ # construct response traits from each linear predictor
    d[paste0("y",i)] = u[i] + unlist(u.list[i]) + unlist(e.list[i])
  }
    
  return (d)

}


level_cor_vcv <- function(vcv, n_resp, part.1 = "animal", part.2 =  "units") {
  
  # index corrmat elements for output
  len_cor <- choose(n_resp,2)
  m <- matrix(1:(n_resp^2),n_resp,n_resp)
  for (i in 1:nrow(m)) {
    for (j in i:(ncol(m))) {
      m[i,j] <- 0
    }
  }
  cols <- unique(as.vector(m));cols <- cols[cols!=0]
  
  # phylogenetic level
  draws <- vcv %>% as_tibble() %>% rename_with(tolower) %>% rename_with(~str_remove_all(., 'trait'))
  cor_draws <- draws %>% dplyr::select(all_of(contains(part.1))) %>% rename_with(~str_remove_all(., paste0('.',part.1)))
  phy_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  par_phy_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  for (i in 1:nrow(cor_draws)){
    cor <- unlist(cor_draws[i,])
    phy_res <- matrix(0, n_resp, n_resp)
    phy_res[] <- cor
    phy_res <- cov2cor(phy_res)
    par_phy_res <- corpcor::cor2pcor(phy_res) # compute partial correlation
    phy_cor[i,] <- as.vector(phy_res)
    par_phy_cor[i,] <- as.vector(par_phy_res)
  }
  phy_cor <- phy_cor %>% as_tibble %>% dplyr::select(all_of(cols)) %>% setNames(names(cor_draws %>% dplyr::select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  par_phy_cor <- par_phy_cor %>% as_tibble %>% dplyr::select(all_of(cols)) %>% setNames(names(cor_draws %>% dplyr::select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  phy_res[lower.tri(phy_res)] <- phy_cor %>% summarise_all(~ mean(.)) %>% slice(1) %>% as.numeric(); phy_res[upper.tri(phy_res)] <- t(phy_res)[upper.tri(phy_res)]; dimnames(phy_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(phy_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  par_phy_res[lower.tri(par_phy_res)] <- par_phy_cor %>% summarise_all(~ mean(.)) %>% slice(1) %>% as.numeric(); par_phy_res[upper.tri(par_phy_res)] <- t(par_phy_res)[upper.tri(par_phy_res)]; dimnames(par_phy_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(par_phy_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  
  # residual level
  draws <- vcv %>% as_tibble() %>% rename_with(tolower) %>% rename_with(~str_remove_all(., 'trait'))
  cor_draws <- draws %>% dplyr::select(all_of(contains(part.2))) %>% rename_with(~str_remove_all(., paste0('.',part.2)))
  ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  par_ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  for (i in 1:nrow(cor_draws)){
    cor <- unlist(cor_draws[i,])
    ind_res <- matrix(0, n_resp, n_resp)
    ind_res[] <- cor
    ind_res <- cov2cor(ind_res)
    par_ind_res <- corpcor::cor2pcor(ind_res) # compute partial correlation
    ind_cor[i,] <- as.vector(ind_res)
    par_ind_cor[i,] <- as.vector(par_ind_res)
  }
  ind_cor <- ind_cor %>% as_tibble %>% dplyr::select(all_of(cols)) %>% setNames(names(cor_draws %>% dplyr::select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  par_ind_cor <- par_ind_cor %>% as_tibble %>% dplyr::select(all_of(cols)) %>% setNames(names(cor_draws %>% dplyr::select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  ind_res[lower.tri(ind_res)] <- ind_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); ind_res[upper.tri(ind_res)] <- t(ind_res)[upper.tri(ind_res)]; dimnames(ind_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(ind_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  par_ind_res[lower.tri(par_ind_res)] <- par_ind_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); par_ind_res[upper.tri(par_ind_res)] <- t(par_ind_res)[upper.tri(par_ind_res)]; dimnames(par_ind_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(par_ind_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))

  result <- list(posteriors = list(phy_cor=phy_cor,par_phy_cor=par_phy_cor,
                                   ind_cor=ind_cor,par_ind_cor=par_ind_cor), 
                 matrices = list(phy_mat = phy_res, par_phy_mat = par_phy_res,
                                 ind_mat = ind_res, par_ind_mat = par_ind_res))
  
  return(result)
  
}

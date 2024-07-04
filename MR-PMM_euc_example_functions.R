# function to calculate phylogenetic signal and phylogenetic and non-phylogenetic correlations 
# between traits from an MR-PMM
mrpmm_cor <- function (fit, y1, y2, phy_group, ind_group, y1_D=0, y2_D=0, summary=F){
  
  phy_cor <- fit$VCV[,paste0("trait",y1,":trait",y2,".",phy_group)] /
    sqrt(fit$VCV[,paste0("trait",y1,":trait",y1,".",phy_group)] * 
           fit$VCV[,paste0("trait",y2,":trait",y2,".",phy_group)])
  ind_cor <- fit$VCV[,paste0("trait",y1,":trait",y2,".",ind_group)] /
    sqrt(fit$VCV[,paste0("trait",y1,":trait",y1,".",ind_group)] * 
           fit$VCV[,paste0("trait",y2,":trait",y2,".",ind_group)])

  y1_sig <- fit$VCV[,paste0("trait",y1,":trait",y1,".",phy_group)] /
                   (fit$VCV[,paste0("trait",y1,":trait",y1,".",phy_group)] + 
                      fit$VCV[,paste0("trait",y1,":trait",y1,".",ind_group)] + y1_D)
  y2_sig <- fit$VCV[,paste0("trait",y2,":trait",y2,".",phy_group)] /
                  (fit$VCV[,paste0("trait",y2,":trait",y2,".",phy_group)] + 
                     fit$VCV[,paste0("trait",y2,":trait",y2,".",ind_group)] + y2_D)
  
  if (summary == T) {
    return(cbind(data.frame(parameter=c(paste0(y1,"_sig"),paste0(y2,"_sig"),"phy_cor","ind_cor")),
                 rbind(round(quantile(y1_sig, probs = c(0.025,0.25,0.5,0.75,0.975)),3),
                       round(quantile(y2_sig, probs = c(0.025,0.25,0.5,0.75,0.975)),3),
                       round(quantile(phy_cor, probs = c(0.025,0.25,0.5,0.75,0.975)),3),
                       round(quantile(ind_cor, probs = c(0.025,0.25,0.5,0.75,0.975)),3))))
  } else {
    
    result <- data.frame(y1_sig,y2_sig,phy_cor,ind_cor)
    names(result) <- c(paste0(y1,"_sig"),paste0(y2,"_sig"),"phy_cor","ind_cor")
    return(result)
  }
}

# function to compute the SE of a weighted mean
# from: https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
wtd.stderror <- function(x, w){
  var <- Hmisc::wtd.var(x, w)
  w <- sum((w / sum(w))^2)
  sqrt(var*w)
}

se <- function(x){sd(x)/sqrt(length(x))}

tree.ultra <- function(tree, data) {
  phy <- keep.tip(tree, data[data$taxon%in%tree$tip.label,]$taxon)
  for (i in 1:3){phy$edge.length[phy$edge.length<=0.00001] <- 0.005
  phy <- multi2di(phy)
  phy <- force.ultrametric(phy, method = "extend", message=F)}
  phy$node.label <- (length(phy$tip.label)+1):(phy$Nnode+length(phy$tip.label))
  return(phy)
}

level_cor_vcv <- function(vcv, n_resp, part.1 = "phylo", part.2 =  "taxon", lambda = NULL) {
  
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
  shrink_par_phy_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  for (i in 1:nrow(cor_draws)){
    cor <- unlist(cor_draws[i,])
    phy_res <- matrix(0, n_resp, n_resp)
    phy_res[] <- cor
    phy_res <- cov2cor(phy_res)
    par_phy_res <- corpcor::cor2pcor(phy_res) # compute partial correlation
    shrink_par_phy_res <- if(is.null(lambda)) corpcor::pcor.shrink(phy_res, verbose = F)
    else(corpcor::pcor.shrink(phy_res, lambda = lambda, verbose = F))# compute shrinkage estimates of partial correlation
    phy_cor[i,] <- as.vector(phy_res)
    par_phy_cor[i,] <- as.vector(par_phy_res)
    shrink_par_phy_cor[i,] <- as.vector(shrink_par_phy_res)
  }
  phy_cor <- phy_cor %>% as_tibble %>% select(all_of(cols)) %>% setNames(names(cor_draws %>% select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  par_phy_cor <- par_phy_cor %>% as_tibble %>% select(all_of(cols)) %>% setNames(names(cor_draws %>% select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  shrink_par_phy_cor <- shrink_par_phy_cor %>% as_tibble %>% select(all_of(cols)) %>% setNames(names(cor_draws %>% select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  phy_res[lower.tri(phy_res)] <- phy_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); phy_res[upper.tri(phy_res)] <- t(phy_res)[upper.tri(phy_res)]; dimnames(phy_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(phy_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  par_phy_res[lower.tri(par_phy_res)] <- par_phy_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); par_phy_res[upper.tri(par_phy_res)] <- t(par_phy_res)[upper.tri(par_phy_res)]; dimnames(par_phy_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(par_phy_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  shrink_par_phy_res[lower.tri(shrink_par_phy_res)] <- shrink_par_phy_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); shrink_par_phy_res[upper.tri(shrink_par_phy_res)] <- t(shrink_par_phy_res)[upper.tri(shrink_par_phy_res)]; dimnames(shrink_par_phy_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(shrink_par_phy_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  
  
  
  # residual level
  draws <- vcv %>% as_tibble() %>% rename_with(tolower) %>% rename_with(~str_remove_all(., 'trait'))
  cor_draws <- draws %>% dplyr::select(all_of(contains(part.2))) %>% rename_with(~str_remove_all(., paste0('.',part.2)))
  ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  par_ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  shrink_par_ind_cor <- matrix(nrow = nrow(cor_draws), ncol = ncol(cor_draws))
  for (i in 1:nrow(cor_draws)){
    cor <- unlist(cor_draws[i,])
    ind_res <- matrix(0, n_resp, n_resp)
    ind_res[] <- cor
    ind_res <- cov2cor(ind_res)
    par_ind_res <- corpcor::cor2pcor(ind_res) # compute partial correlation
    shrink_par_ind_res <- if(is.null(lambda)) corpcor::pcor.shrink(ind_res, verbose = F)
    else(corpcor::pcor.shrink(ind_res, lambda = lambda, verbose = F)) # compute shrinkage estimates of partial correlation
    ind_cor[i,] <- as.vector(ind_res)
    par_ind_cor[i,] <- as.vector(par_ind_res)
    shrink_par_ind_cor[i,] <- as.vector(shrink_par_ind_res)
  }
  ind_cor <- ind_cor %>% as_tibble %>% select(all_of(cols)) %>% setNames(names(cor_draws %>% select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  par_ind_cor <- par_ind_cor %>% as_tibble %>% select(all_of(cols)) %>% setNames(names(cor_draws %>% select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  shrink_par_ind_cor <- shrink_par_ind_cor %>% as_tibble %>% select(all_of(cols)) %>% setNames(names(cor_draws %>% select(all_of(cols)))) # %>% mcmc_intervals(prob_outer = 0.95)
  ind_res[lower.tri(ind_res)] <- ind_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); ind_res[upper.tri(ind_res)] <- t(ind_res)[upper.tri(ind_res)]; dimnames(ind_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(ind_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  par_ind_res[lower.tri(par_ind_res)] <- par_ind_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); par_ind_res[upper.tri(par_ind_res)] <- t(par_ind_res)[upper.tri(par_ind_res)]; dimnames(par_ind_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(par_ind_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  shrink_par_ind_res[lower.tri(shrink_par_ind_res)] <- shrink_par_ind_cor %>% summarise(across(everything(), ~ mean(.))) %>% slice(1) %>% as.numeric(); shrink_par_ind_res[upper.tri(shrink_par_ind_res)] <- t(shrink_par_ind_res)[upper.tri(shrink_par_ind_res)]; dimnames(shrink_par_ind_res)[[1]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])));dimnames(shrink_par_ind_res)[[2]] <- gsub(":","",gsub("trait","",gsub("[^:]+$","",colnames(vcv)[1:n_resp])))
  
  
  result <- list(posteriors = list(phy_cor=phy_cor,par_phy_cor=par_phy_cor,shrink_par_phy_cor=shrink_par_phy_cor,
                                   ind_cor=ind_cor,par_ind_cor=par_ind_cor,shrink_par_ind_cor=shrink_par_ind_cor), 
                 matrices = list(phy_mat = phy_res, par_phy_mat = par_phy_res, shrink_par_phy_mat = shrink_par_phy_res,
                                 ind_mat = ind_res, par_ind_mat = par_ind_res, shrink_par_ind_mat = shrink_par_ind_res))
  
  return(result)
  
}

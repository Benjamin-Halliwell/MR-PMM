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

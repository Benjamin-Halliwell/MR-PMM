library(ips);library(ggpubr);library(tidyr);library(data.table);library(caper);library(RRphylo)
library(stringr);library(stringdist);library(phangorn);library(corrplot);library(grid)
library(RColorBrewer);library(phytools);library(dplyr);library(plyr);library(rr2);library(motmot)
library(phylolm);library(MCMCglmm);library(diversitree);library(ggtree);library(devtools)
library(BiocManager);library(tidytree);library(cowplot);library(ggpubr);library(patchwork)
library(surface);library(ouch);library(geiger);library(igraph);library(Rphylopars);library(rlist)
require(dplyr);library(brms);library(future)


force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}


### DATA AND TREES ####

# d.master <- read.csv("d.master.csv")
# d.master <- d.master[d.master$num_occ>2,]
d.master <- read.csv("d.master_V4.csv")


## TREES
# tree <- read.tree(file = 'Eucalypts_ML1_dated_r8s.tre')
tree.ML1 <- read.tree(file="tree.ML1.tre")
tree.ML2 <- read.tree(file="tree.ML2.tre")
tree.Bayes <- read.tree(file="tree.Bayes.tre")

## Sequence data from Thornhill
# seq <- read.fas("Eucalypts_final_PD.fas")

# resolve polytomies and make clock-like
tree.ML1 <- multi2di(tree.ML1, random = F)
tree.ML2 <- multi2di(tree.ML2, random = F)
tree.Bayes <- multi2di(tree.Bayes, random = F)
tree.ML1$node.label <- (length(tree.ML1$tip.label)+1):(tree.ML1$Nnode+length(tree.ML1$tip.label))
tree.ML2$node.label <- (length(tree.ML1$tip.label)+1):(tree.ML1$Nnode+length(tree.ML1$tip.label))
tree.Bayes$node.label <- (length(tree.ML1$tip.label)+1):(tree.ML1$Nnode+length(tree.ML1$tip.label))
tree.ML1.u <- force.ultrametric(tree.ML1)
tree.ML2.u <- force.ultrametric(tree.ML2)
tree.Bayes.u <- force.ultrametric(tree.Bayes)

# subset data to Bayes phy
d.bayes <- d.master[d.master$taxon %in% tree.Bayes$tip.label,]

### END ####


### PMM ####

# trees
phy1 <- tree.ML1.u
phy2 <- tree.ML2.u
phy3 <- tree.Bayes.u

# set plot view
par(mar=(c(2,2,2,2)))

# resolve zero branch lengths
for (i in 1:3){phy1$edge.length[phy1$edge.length<=0.00001] <- 0.0001
phy1 <- force.ultrametric(phy1)}

for (i in 1:3){phy2$edge.length[phy2$edge.length<=0.00001] <- 0.0001
phy2 <- force.ultrametric(phy2)}

for (i in 1:3){phy3$edge.length[phy3$edge.length<=0.00001] <- 0.0001
phy3 <- force.ultrametric(phy3)}


## SUBSET TO SYMPHYOMYRTUS AND EUCALYPTUS ##
d.subgen <- d.master[d.master$subgenus=="Symphyomyrtus"|d.master$subgenus=="Eucalyptus",]
row.names(d.subgen) <-  1:nrow(d.subgen)
d.subgen$subgenus <- as.factor(d.subgen$subgenus)


tree.subgen <- keep.tip(phy1, d.subgen$taxon)
for (i in 1:3){tree.subgen$edge.length[tree.subgen$edge.length<=0.00001] <- 0.005
tree.subgen <- force.ultrametric(tree.subgen)}



## COMPLETE CASES FOR SE
# d.sub <- d.subgen[!is.na(d.subgen$SLA.log)&                   
#                     !is.na(d.subgen$leaf_N_per_dry_mass)&
#                     !is.na(d.subgen$leaf_delta13C)&
#                     !is.na(d.subgen$se.specific_leaf_area)&
#                     !is.na(d.subgen$se.leaf_N_per_dry_mass)&
#                     !is.na(d.subgen$se.leaf_delta13C),]

## COMPLETE CASES FOR TRAITS
d.sub <- d.subgen[!is.na(d.subgen$SLA)&
                    !is.na(d.subgen$leaf_N_per_dry_mass)&
                    !is.na(d.subgen$leaf_delta13C),]

# subset data and other trees to Bayes tree
d.sub <- d.sub[d.sub$taxon %in% phy3$tip.label,]
phy1 <- keep.tip(phy1, d.sub$taxon)
phy2 <- keep.tip(phy2, d.sub$taxon)
phy3 <- keep.tip(phy3, d.sub$taxon)
names(d.sub)[names(d.sub)=="se.specific_leaf_area"] <- "se.SLA"
d.sub$obs <- 1:nrow(d.sub)

# ensure all SE are positive by adding min SE observed for trait to any that equal 0
d.sub[!is.na(d.sub$se.SLA) & d.sub$se.SLA==0,]$se.SLA <- min(d.sub[!is.na(d.sub$se.SLA) & d.sub$se.SLA!=0,]$se.SLA)
d.sub[!is.na(d.sub$se.leaf_N_per_dry_mass) & d.sub$se.leaf_N_per_dry_mass==0,]$se.leaf_N_per_dry_mass <- min(d.sub[!is.na(d.sub$se.leaf_N_per_dry_mass) & d.sub$se.leaf_N_per_dry_mass!=0,]$se.leaf_N_per_dry_mass)
d.sub[!is.na(d.sub$se.leaf_delta13C) & d.sub$se.leaf_delta13C==0,]$se.leaf_delta13C <- min(d.sub[!is.na(d.sub$se.leaf_delta13C) & d.sub$se.leaf_delta13C!=0,]$se.leaf_delta13C)

# for any NAs, make SE the 90th pctl. of SE for that trait
d.sub[is.na(d.sub$se.SLA),]$se.SLA <- quantile(d.sub[!is.na(d.sub$se.SLA),]$se.SLA, 0.9)
d.sub[is.na(d.sub$se.leaf_N_per_dry_mass),]$se.leaf_N_per_dry_mass <- quantile(d.sub[!is.na(d.sub$se.leaf_N_per_dry_mass),]$se.leaf_N_per_dry_mass, 0.9)
d.sub[is.na(d.sub$se.leaf_delta13C),]$se.leaf_delta13C <- quantile(d.sub[!is.na(d.sub$se.leaf_delta13C),]$se.leaf_delta13C, 0.9)


# log-transform SLA and compute first-order standard error estimates (i.e. relative error) via delta method (se.SLA/SLA)
d.sub <- d.sub %>% as_tibble %>% mutate(log_SLA = log(SLA), se_log_SLA = se.SLA/SLA,
                                        N = leaf_N_per_dry_mass, se_N = se.leaf_N_per_dry_mass,
                                        d13C = leaf_delta13C, se_d13C = se.leaf_delta13C)
# subset to relevant columns
d.sub2 <- d.sub %>% select(animal, log_SLA, se_log_SLA,N, se_N,d13C,se_d13C)


# ## SAVED FITS
# b.0 <- readRDS("b.0.rds")
# b.1 <- readRDS("b.1.rds")
# b.2 <- readRDS("b.2.rds")
# b.3 <- readRDS("b.3.rds")

## brms
A.mat <- ape::vcv.phylo(phy1, corr = T)

## MOD 0 - raw trait scale
bf_y1 <- bf(SLA | resp_se(se.SLA, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A)))
bf_y2 <- bf(leaf_N_per_dry_mass | resp_se(se.leaf_N_per_dry_mass, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A)))
bf_y3 <- bf(leaf_delta13C | resp_se(se.leaf_delta13C, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A)))

b.0 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
           data = d.sub, 
           family = gaussian(), 
           data2 = list(A = A.mat),
           control=list(adapt_delta = 0.85, max_treedepth = 12),
           cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.0, file = "b.0.rds")
# b.0
# plot(b.0) # posteriors for phy cor terms skewed and stuck around 0 or 1

# pp checks show SLA needs to be transform to avoid boundary issues
# p1 <- pp_check(b.0, resp = "SLA", ndraws = 100) + ggtitle("SLA") # prediction poor for SLA
# p2 <- pp_check(b.0, resp = "leafNperdrymass", ndraws = 100) + ggtitle("N")
# p3 <- pp_check(b.0, resp = "leafdelta13C", ndraws = 100) + ggtitle("d13C")
# ggarrange(p1,p2,p3, nrow = 1, ncol = 3)


# fit with log transformed SLA instead

## phy1
bf_y1 <- bf(log_SLA | resp_se(se_log_SLA, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A))) 
bf_y2 <- bf(N | resp_se(se_N, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A))) 
bf_y3 <- bf(d13C | resp_se(se_d13C, sigma = TRUE) ~ 1 + (1|b|gr(animal, cov = A))) 

b.1 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
            data = d.sub2, 
            family = gaussian(), 
            data2 = list(A = A.mat),
            control=list(adapt_delta = 0.85, max_treedepth = 12),
            cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.1, file = "b.1.rds")
# readRDS(file = "b.1.rds")

# # pp checks now look good. COntinue fit for other trees
# p1 <- pp_check(b.1, resp = "logSLA", ndraws = 100) + ggtitle("SLA") # prediction poor for SLA
# p2 <- pp_check(b.1, resp = "N", ndraws = 100) + ggtitle("N")
# p3 <- pp_check(b.1, resp = "d13C", ndraws = 100) + ggtitle("d13C")
# ggarrange(p1,p2,p3, nrow = 1, ncol = 3)


## phy2
A.mat2 <- ape::vcv.phylo(phy2, corr = T)
b.2 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
             data = d.sub2, 
             family = gaussian(), 
             data2 = list(A = A.mat2),
             control=list(adapt_delta = 0.85, max_treedepth = 12),
             cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.2, file = "b.2.rds")


## phy3
A.mat3 <- ape::vcv.phylo(phy3, corr = T)
b.3 <- brm(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
             data = d.sub2, 
             family = gaussian(), 
             data2 = list(A = A.mat3),
             control=list(adapt_delta = 0.85, max_treedepth = 12),
             cores = 4, chains = 4, iter = 6000, thin = 3)
saveRDS(b.3, file = "b.3.rds")


## COMBINE MODELS ##

## results similar across 3 available trees so combine
# b.mod <- combine_models(b.1, b.2, b.3)
# saveRDS(b.mod, file = "b.mod.rds")
b.mod <- readRDS(file = "b.mod.rds")
b.mod

# worth fitting a skew normal to SLA?
p1 <- pp_check(b.mod, resp = "logSLA", ndraws = 100) + ggtitle("SLA") # prediction poor for SLA
p2 <- pp_check(b.mod, resp = "N", ndraws = 100) + ggtitle("N")
p3 <- pp_check(b.mod, resp = "d13C", ndraws = 100) + ggtitle("d13C")
ggarrange(p1,p2,p3, nrow = 1, ncol = 3)




## BRMS Multiple ##

# make list of data and trees
t.list <- list(tree.sub1,tree.sub2)
t.list <- lapply(t.list, ape::vcv.phylo)
t.list <- list(list(A=t.list[[1]]),list(A=t.list[[2]]))
d.list <- list(d.sub, d.sub)

# fit over trees and parallel process
plan(multiprocess)
brm_multiple(bf_y1 + bf_y2 + bf_y3 + set_rescor(TRUE),
             data = d.list,
             data2 = t.list,
             family = gaussian(), 
             control=list(adapt_delta = 0.85, max_treedepth = 12),
             cores = 4, chains = 4, iter = 6000, thin = 3)


# ggplot(d.sub, aes(x=leaf_N_per_dry_mass,y=leaf_delta13C,col=subgenus)) + geom_point() + theme_classic() + geom_smooth(method = lm)


### END ####





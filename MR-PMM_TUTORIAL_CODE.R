library("MCMCglmm");library("brms");library("tidyverse");library("ape");
library("phytools");library("MASS");library("plyr");library("bindata");library('phangorn')


##---------------------------- ~ --------------------------------##

#### DEFINE SIMULATION PARAMETERS ####

# simulate tree
t.toy <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=250, extinct=FALSE) # may be best not to use pure birth trees (Adams and Collyer 2018)
t.toy <- multi2di(t.toy)
plot(t.toy)
# save(t.toy, file="t.toy")

# Create VCV matrix from tree, scale to correlation matrix for BRMS
A.mat <- ape::vcv.phylo(t.toy, corr = T) # scale with corr = T

# ## alternative method 
# inv.phylo <- MCMCglmm::inverseA(t.toy, nodes = "TIPS", scale = TRUE)
# A <- solve(inv.phylo$Ainv)
# rownames(A) <- rownames(inv.phylo$Ainv)
# A.mat <- A

# number of traits
k = 2

# intercepts for each trait - plogis(beta[1]) to get on probability scale
beta = c(1,2)

## B is the phylogenetic trait-level VCV matrix (Sigma_{phy} in the HTML). B specifies the phylogenetic
## variance in each trait (diagonals) as well as the phylogenetic covariance between traits (off-diagnonals).
## Each unique element in B is to be estimated by the model.
sig.B <- c(b11 = 0.4, b22 = 0.6) # standard deviation of diagonal components (i.e. sqrt of the phylogenetic variance for each trait)
b12_rho = 0.5

Bcor <- matrix(c(c(1,b12_rho), # correlation matrix
                 c(b12_rho,1)),k,k, byrow = T) 

B <- matrix(kronecker(sig.B, sig.B),k,k)*Bcor # VCV (point-wise product). 
B
# N.B. Kronecker used here just for matrix formatting. Do not confuse with kronecker discussed in the text / tutorial.

## C is the residual trait-level VCV matrix (Sigma_{res} in the HTML).
sig.C <- c(c11 = 0.1, c22 = 0.1) # standard deviation of diagonal components (i.e. sqrt of the residual variance for each trait)
c12_rho = 0 # off-diagonal correlation coefficients

Ccor <- matrix(c(c(1,c12_rho), # correlation matrix
                 c(c12_rho,1)),k,k, byrow = T) 

C <- matrix(kronecker(sig.C, sig.C),k,k)*Ccor # VCV
C

# The phylogenetic taxon-level VCV matrix, A, is supplied. It is therefore treated as fixed and known without error, although may represent a transformation (see above).
A = A.mat

# number of species
n = nrow(A); n

# identity matrix
I = diag(n)

# simulate phylogenetic random effects for all traits as one draw from a MVN:
a <- mvrnorm(1,rep(0,n*k),kronecker(B,A)) # In the Kronecker, trait-level covariance captured by B, taxon-level covariance captured by A. 

# extract random effects for each trait from the vector, a.
a1 <- a[1:n]
a2 <- a[1:n + n]

# simulate residuals (on the link scale)
e <- mvrnorm(1,rep(0,n*k),kronecker(C,I))
e1 <- e[1:n]
e2 <- e[1:n + n]

# construct response traits from each linear predictor
g1 <- beta[1] + a1 + e1 # gaussian
g2 <- beta[2] + a2 + e2
b1 <- rbinom(n,1,plogis(beta[1] + a1 + e1)) # binomial
b2 <- rbinom(n,1,plogis(beta[2] + a2 + e2))

# generate df
species <- t.toy$tip.label
d.toy <- data.frame(species,g1,g2,b1,b2)
d.toy$animal <- d.toy$species # "animal" is a reserved term in MCMCglmm used to identify taxa in quantitative and phylogenetic models
d.toy$obs <- 1:nrow(d.toy)
head(d.toy)

## plot gaussians
# plot(d.toy$g1, d.toy$g2)
# table binomials
table(d.toy$b1, d.toy$b2)
table(d.toy$b1, d.toy$b2)/250 # prop

## take a look at probabilities on the inverse-link scale
# plogis(beta[1])
# plogis(beta[1] + a1 + e1)[1:10]

### END ####




##-------------------- UNIVARIATE GAUSSIAN ----------------------##


## MCMCglmm 
p <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
          R = list(V = 1, nu = 0.002))

m.1<-MCMCglmm(g1 ~ 1,
              random = ~animal,
              pedigree=t.toy, 
              family = c("gaussian"), 
              data = d.toy,  prior = p,
              nitt=600000, burnin=100000, thin=500,
              pr=TRUE,verbose = FALSE)
summary(m.1)

## DIAGNOSTICS
par(mar=(c(2,2,2,2)))
plot(m.1$VCV)
# calculate autocorrelation among samples. All values above lag 0 should be close to 0
autocorr(m.1$VCV)


## BRMS
b.1 <- brm(
  g1 ~ 1 + (1|gr(animal, cov = A)), 
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 2000, thin = 2
)
b.1

## DIAGNOSTICS
b.1 %>% plot
b.1 %>% pp_check(nsamples = 100)
pairs(b.1) # sd_animal and sigma are trading off


## COMPARISON BETWEEN METHODS

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

# N.B. Intercept estimates are directly comparable between MCMCglmm and BRMS but the variance components reported
# by MCMCglmm need to be sqrt() to compare to the SDs reported by BRMS (as well as the generating parameters for our sims).

## Parameter estimates:
# intercept
summary(m.1)$solutions
summary(b.1)[["fixed"]]
# phylogenetic variance (sig.B)
sqrt(summary(m.1)$Gcovariances) # N.B. effective sample size is now incorrect
# sd(m.1$VCV[, 1]) # or take sd directly?
summary(b.1)[["random"]]
# residual variance (sig.C)
sqrt(summary(m.1)$Rcovariances)
summary(b.1)[["spec_pars"]]



##------------------- MULTIVARIATE GAUSSIAN ---------------------##

m.1 <- readRDS("m.1.rds")
b.1 <- readRDS("b.1.rds")

### FIT MODELS ####

### MCMCglmm
p <- list(R = list(R1=list(V = diag(2), nu=2)),
          G = list(G1=list(V = diag(2), nu=2)))
m.1 <- MCMCglmm(cbind(g1, g2) ~ trait-1, 
                random = ~us(trait):animal, 
                rcov = ~us(trait):units, 
                pedigree=t.toy,
                family = c("gaussian","gaussian"), 
                nodes="ALL", data = d.toy, prior=p,
                nitt=600000, burnin=100000, thin=500, 
                pr=TRUE,verbose = FALSE)

## DIAGNOSTICS
# view chain
# plot(m.1$VCV)
# calculate autocorrelation among samples
# autocorr(m.1$VCV)


### BRMS
b.1 <- brm(
  mvbind(g1, g2) ~ (1|p|gr(animal, cov = A)), 
  data = d.toy,
  data2 = list(A = A.mat),
  family = gaussian(),
  cores = 4,
  chains = 4, iter = 2000, thin = 2
)
b.1
## DIAGNOSTICS
pairs(b.1)
b.1 %>% plot
b.1 %>% pp_check(resp = "g1",ndraws = 100)
b.1 %>% pp_check(resp = "g2",ndraws = 100)
conditional_effects(b.1)


## COMPARISON BETWEEN METHODS
# N.B. Location estimates are directly comparable between MCMCglmm and BRMS (assuming no variances have fixed OR fixed at same level). However, for variance 
# components, MCMCglmm reports (co)variances while brms reports standard deviations and trait-level correlations. Therefore, some re-scaling is necessary to 
# to compare results between model fits. We have chosen to re-scale the estimates from MCMCglmm, as stdevs and corrs are more natural to interpret.

# summary outputs
summary(m.1) # MCMCglmm
summary(b.1) # brms

# intercepts
summary(m.1)$solutions
summary(b.1)[["fixed"]]

# phylogenetic variances (elements of sig.B) and phylogenetic correlation (b12_rho)  # Q. SQRT OF CI CORRECT CI FOR VARIANCES? NO, write function to do this and make into df
sqrt(summary(m.1)$Gcovariances)[c(1,nrow(summary(m.1)$Gcovariances)),];data.frame(row.names = "traitg1:traitg2.animal", post.mean = mean(m.1$VCV[,"traitg1:traitg2.animal"]/sqrt(m.1$VCV[,"traitg1:traitg1.animal"]*m.1$VCV[,"traitg2:traitg2.animal"])))
summary(b.1)[["random"]]

# residual covariances (elements of sig.C) and residual correlation (c12_rho) 
sqrt(summary(m.1)$Rcovariances[c(1,nrow(summary(m.1)$Gcovariances)),]);data.frame(row.names = "traitg1:traitg2.units", post.mean = mean(m.1$VCV[,"traitg1:traitg2.units"]/sqrt(m.1$VCV[,"traitg1:traitg1.units"]*m.1$VCV[,"traitg2:traitg2.units"])))
summary(b.1)[["spec_pars"]];summary(b.1)[["rescor_pars"]]


### END ####

##------------- MULTIVARIATE (GAUSSIAN, BERNOULLI) --------------##

m.2a <- readRDS("m.2a.rds")
m.2b <- readRDS("m.2b.rds")
b.2 <- readRDS("b.2.rds")

### FIT MODELS ####

### MCMCglmm

# Poor mixing of random effects for binary traits with default priors. Use parameter expanded priors to improve mixing.
# Fix residual variance of binary trait at two different levels (1 and 5) to check if this further improves mixing.
# If so, estimates can be rescaled to a specified level of variance for comparison (see below).

p2a <- list(B = list(mu=c(0,0), V=diag(c(1,1+pi^2/3))),
            R = list(R1=list(V = diag(2), nu=2, fix=2)),
            G = list(G1=list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V = diag(c(1,1000)))))

p2b <- list(B = list(mu=c(0,0), V=diag(c(1,1+pi^2/3))),
            R = list(R1=list(V = diag(c(1,5)), nu=2, fix=2)),
            G = list(G1=list(V = diag(2), nu=3, alpha.mu = c(0,0), alpha.V = diag(c(1,1000)))))

m.2a <- MCMCglmm(cbind(g1, b2) ~ trait-1, 
                random = ~us(trait):animal, 
                rcov = ~us(trait):units, # using ~us() while fixing residual variance for b2 and residual covariance between g1 and b2 in the prior. Could explore corg() to fix var but still estimate res cov on link scale?
                pedigree=t.toy,
                family = c("gaussian","categorical"), 
                nodes="ALL", data = d.toy, prior=p2a,
                nitt=600000, burnin=100000, thin=500, 
                pr=TRUE,verbose = FALSE) 
saveRDS(m.2a, "m.2a.rds")

m.2b <- MCMCglmm(cbind(g1, b2) ~ trait-1, random = ~us(trait):animal, rcov = ~us(trait):units, pedigree=t.toy,family = c("gaussian","categorical"), nodes="ALL", data = d.toy, prior=p2b,nitt=600000, burnin=100000, thin=500, pr=TRUE,verbose = FALSE) 
saveRDS(m.2b, "m.2b.rds")


## DIAGNOSTICS
# view chain
plot(m.2a$VCV)
plot(m.2b$VCV)

# calculate autocorrelation among samples
autocorr(m.2a$VCV)
autocorr(m.2b$VCV)

summary(m.2a)
summary(m.2b)


#----------------------------------------------------------#


### BRMS

bf_g1 <- bf(g1 ~ 1 + (1|a|gr(animal, cov = A)) + (1|obs), sigma = 0.001) + gaussian() # fix residual error at very low level to force residual error into 1|obs
bf_b2 <- bf(b2 ~ 1 + (1|a|gr(animal, cov = A)) + (1|obs)) + bernoulli() # (1|obs) rather than (1|b|obs) so residual covariance not estimated

b.2 <- brm(
  bf_g1 + bf_b2 + set_rescor(FALSE), 
  data = d.toy,
  data2 = list(A = A.mat),
  prior = set_prior("constant(1)", class = "sd", group = "obs", resp = "b2"),
  family = gaussian(), # THIS OK TO STAY GAUSSIAN?
  cores = 4,
  chains = 4, iter = 2000, thin = 2
)
saveRDS(b.2, "b.2.rds")

## DIAGNOSTICS
b.2 %>% plot
b.2 %>% pp_check(resp = "g1",ndraws = 100)
b.2 %>% pp_check(resp = "g2",ndraws = 100)
conditional_effects(b.2)

## COMPARE

# summary outputs
summary(m.2)
summary(b.2)

# intercept
summary(m.2)$solutions
summary(b.2)[["fixed"]]

# phylogenetic variances (elements of sig.B) and phylogenetic correlation (b12_rho)  # Q. SQRT OF CI CORRECT CI FOR VARIANCES?
sqrt(summary(m.1)$Gcovariances)[c(1,nrow(summary(m.1)$Gcovariances)),];data.frame(row.names = "traitg1:traitb2.animal", post.mean = mean(m.1$VCV[,"traitg1:traitb2.animal"]/sqrt(m.1$VCV[,"traitg1:traitg1.animal"]*m.1$VCV[,"traitb2:traitb2.animal"])))
summary(b.1)[["random"]]

# residual covariances (elements of sig.C) and residual correlation (c12_rho) 
sqrt(summary(m.1)$Rcovariances[c(1,nrow(summary(m.1)$Gcovariances)),]);data.frame(row.names = "traitg1:traitb2.units", post.mean = mean(m.1$VCV[,"traitg1:traitb2.units"]/sqrt(m.1$VCV[,"traitg1:traitg1.units"]*m.1$VCV[,"traitb2:traitb2.units"])))
summary(b.1)[["spec_pars"]];summary(b.1)[["rescor_pars"]]


### END ####

##------------- MULTIVARIATE (BERNOULLI, BERNOULLI) -------------##

### FIT MODELS ####

### MCMCglmm 

# For binary responses, the residual variance is not identifiable. Therefore, for family = categorical residual variances need to be fixed
# in the prior specification (ref: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q3/023966.html) and rcov = ~idh() specified to signify that residual 
# covariance should not be estimated. Another option is to use rcov = ~corg() with family = threshold to fix residual variances to 1 while still allowing estimation 
# of residual covariance.

p3a <- list(G = list(G1=list(V=diag(2), nu=1000, alpha.mu=c(0,0), alpha.V=diag(c(1,1)))), # Chi^2 - parameter expanded G priors increase ESS of intercept and phylogenetic variance
            R = list(V=diag(2), nu=0, fix=1))


# p3a <- list(G = list(G1=list(V=diag(2), nu=3, alpha.mu=c(0,0), alpha.V=diag(c(1000,1000)))), # scaled Fisher - parameter expanded G priors increase ESS of intercept and phylogenetic variance
#             R = list(V=diag(2), nu=0, fix=1))

p3b <- list(G = list(G1=list(V=diag(c(5,5)), nu=3, alpha.mu=c(0,0), alpha.V=diag(c(1000,1000)))), # scaled Fisher - parameter expanded G priors increase ESS of intercept and phylogenetic variance
            R = list(V=diag(2), nu=0, fix=1))


m.3a<-MCMCglmm(cbind(b1, b2) ~ trait-1,
              random = ~us(trait):animal, 
              rcov = ~us(trait):units,
              family = c("categorical","categorical"),
              pedigree = t.toy,
              data = d.toy, 
              prior = p3a,
              nitt=600000, burnin=100000, thin=500, 
              pr = T, verbose = F)
summary(m.3a)

m.3b<-MCMCglmm(cbind(b1, b2) ~ trait-1,
               random = ~us(trait):animal, 
               rcov = ~us(trait):units,
               family = c("categorical","categorical"),
               pedigree = t.toy,
               data = d.toy, 
               prior = p3b,
               nitt=600000, burnin=100000, thin=500, 
               pr = T, verbose = F)


## THRESHOLD MODEL 
# For bivariate problems, Hadfield recommends use of the 'threshold' family with probit link and a scaled Fisher prior
# on the random effects. Refs:
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q3/023966.html
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q1/021875.html

p3c <- list(G = list(G1=list(V=diag(c(5,5)), nu=3, alpha.mu=c(0,0), alpha.V=diag(c(1000,1000)))),
            R = list(V=diag(2), nu=0))

m.3c<-MCMCglmm(cbind(b1, b2) ~ trait-1,
              random = ~us(trait):animal, 
              rcov = ~corg(trait):units,
              family = c("threshold","threshold"),
              pedigree = t.toy,
              data = d.toy, 
              prior = p3c,
              nitt=600000, burnin=100000, thin=500,
              pr = T, verbose = F)

## DIAGNOSTICS
plot(m.3a$VCV)
plot(m.3b$VCV)
plot(m.3c$VCV)

autocorr(m.3a$VCV)
autocorr(m.3b$VCV)
autocorr(m.3c$VCV)


#---------------------------------------------------------------------#

### BRMS

# The residual variance of a binomial variable is non-identifiable. Therefore, when attempting to estimate
# group level (co)variances in a MV-PMM involving a binary trait, we need to fix the residual variance in order to... 
# Unlike MCMCglmm, residual variance is not fit by default for a binomial response in brms. We must specify it manually:
#   + (1|b|obs) models residual (co)variance between traits
#   + (1|obs) models only residual variance (additive overdispersion in the binomial case)
# By specifying + (1|obs) AND fixing the residual variance of both responses via the prior specification in set_prior(),
# we are able to replicate the MCMCglmm model with fixed residual (co)variances.

# see which priors brms is using by default for our desired model
priors_b.1<-get_prior(mvbind(b1,b2) ~ 1 + (1|a|gr(animal, cov = A)) + (1|obs),
                      data = d.toy,
                      data2 = list(A = A.mat),
                      family = bernoulli())  
priors_b.1


b.3 <- brm(
  bf(mvbind(b1,b2) ~ 1 + (1|a|gr(animal, cov = A)) + (1|obs)), # + bernoulli() redundant if family = bernoulli(), however must specify separate bf() if prob dist of responses differs
  data = d.toy,
  prior = set_prior("constant(1)", class = "sd", group = "obs", resp = c("b1","b2")),
  family = bernoulli(),
  control = list(adapt_delta = 0.95), # adapt_delta increased to reduce divergent transitions
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 2000, thin = 2
)
saveRDS(b.3, "b.3.rds")

## DIAGNOSTICS
b.3 %>% plot
b.3 %>% pp_check(nsamples = 100)
pairs(b.3)

## COMPARISON
# summary outputs
summary(m.3a)
summary(b.3)
# intercept
summary(m.3a)$solutions
summary(b.3)[["fixed"]]
# phylogenetic (co)variance
sqrt(abs(summary(m.3a)$Gcovariances)) # N.B. remember sign of covariance
summary(b.3)[["random"]]


### END ####

##------ MULTIVARIATE GAUSSIAN, MULTIPLE GROUPING FACTORS -------##

### FIT MODELS ####

# For some analyses, we may wish to model phylogenetic and residual (co)variances separately for responses classified into different levels of a grouping factor
# e.g. estimate the phylogenetic correlation between seed size and asymptotic height separately for woody and herbaceous plant species

# create grouping factor
d.toy$group <- sample(c("L1", "L2"), size = nrow(d.toy), replace = T)

## MCMCglmm
p4=list(G = list(G1=list(V = diag(2), nu = 1.002),
                 G2=list(V = diag(2), nu = 1.002)),
        R = list(R1=list(V = diag(2), nu = 1.002),
                 R2=list(V = diag(2), nu = 1.002)))

m.4<-MCMCglmm(cbind(g1, g2) ~ trait:group-1, # fit separate intercepts for each trait for each level of group
              random = ~us(at.level(group,'L1'):trait):animal + us(at.level(group,'L2'):trait):animal, # covariance matrices specified for each level of group
              rcov = ~us(at.level(group,'L1'):trait):units + us(at.level(group,'L2'):trait):units,
              family = c("gaussian","gaussian"), 
              pedigree=t.toy, data = d.toy, prior=p4, 
              nitt=600000, burnin=100000, thin=500,
              pr=TRUE,verbose = FALSE)
summary(m.4)

## DIAGNOSTICS
plot(m.4$VCV)
autocorr(m.4$VCV)

#--------------------------------------------------------#

## BRMS
b.4 <- brm(
  mvbind(g1, g2) ~ (1|p|gr(by=group, animal, cov = A)), # grouping factors easier to specify in brms: 'by=' compared with 'at.level()' coding in MCMCglmm
  data = d.toy, 
  family = gaussian(), 
  data2 = list(A = A.mat),
  cores=4,
  chains=4, iter = 6000, thin = 3
)
summary(b.4)
b.4 %>% plot
b.4 %>% pp_check(resp = "g1",nsamples = 100)
b.4 %>% pp_check(resp = "g2",nsamples = 100)


## COMPARISON

# summary outputs
summary(m.4)
summary(b.4)

# location effects
summary(m.4)$solutions
summary(b.4)[["fixed"]]

# phylogenetic variances (elements of sig.B)
sqrt(summary(m.4)$Gcovariances)[c(1,nrow(summary(m.4)$Gcovariances)),]
summary(b.4)[["random"]]

## phylogenetic correlation (b12_rho)
# L1
posterior.mode(m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg2.animal']/sqrt((m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg1.animal']*m.4$VCV[,'at.level(group, "L1"):traitg2:at.level(group, "L1"):traitg2.animal'])))
round(HPDinterval(m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg2.animal']/sqrt((m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg1.animal']*m.4$VCV[,'at.level(group, "L1"):traitg2:at.level(group, "L1"):traitg2.animal']))),2)[1:2]
# L2
posterior.mode(m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg2.animal']/sqrt((m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg1.animal']*m.4$VCV[,'at.level(group, "L2"):traitg2:at.level(group, "L2"):traitg2.animal'])))
round(HPDinterval(m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg2.animal']/sqrt((m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg1.animal']*m.4$VCV[,'at.level(group, "L2"):traitg2:at.level(group, "L2"):traitg2.animal']))))

summary(b.4)[["random"]]

# residual covariances (elements of sig.C)
sqrt(summary(m.4)$Rcovariances[c(1,nrow(summary(m.4)$Gcovariances)),])
summary(b.4)[["spec_pars"]]

## residual correlation (c12_rho)
# L1
posterior.mode(m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg2.units']/sqrt((m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg1.units']*m.4$VCV[,'at.level(group, "L1"):traitg2:at.level(group, "L1"):traitg2.units'])))
round(HPDinterval(m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg2.units']/sqrt((m.4$VCV[,'at.level(group, "L1"):traitg1:at.level(group, "L1"):traitg1.units']*m.4$VCV[,'at.level(group, "L1"):traitg2:at.level(group, "L1"):traitg2.units']))),2)[1:2]
# L2
posterior.mode(m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg2.units']/sqrt((m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg1.units']*m.4$VCV[,'at.level(group, "L2"):traitg2:at.level(group, "L2"):traitg2.units'])))
round(HPDinterval(m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg2.units']/sqrt((m.4$VCV[,'at.level(group, "L2"):traitg1:at.level(group, "L2"):traitg1.units']*m.4$VCV[,'at.level(group, "L2"):traitg2:at.level(group, "L2"):traitg2.units']))))
summary(b.4)[["rescor_pars"]]
                                         
### END ####

##---------------------------- ~ --------------------------------##

### TO ADD ####

## INTRACLASS CORRELATION COEFFICIENT (PHYLO SIGNAL FOR BINARY TRAITS) ##

# For binary traits, group level variance only makes sense in context of the residual variance.
# Therefore phylogenetic heritability is not an appropriate measure of signal.
# Instead, use intra-class correlation coefficient, rescaling the phylogenetic variance 
# according to the assumed residual variance (fixed at 1)
head(m.3a$VCV)

IC.b1 <- m.3a$VCV[, 1]/(m.3a$VCV[, 1] + 1 + pi^2/3) # 1 = fixed additive overdispersion, pi^2/3 = family specific variance
mean(IC.b1);HPDinterval(IC.b1)
IC.b2 <- m.3a$VCV[, 4]/(m.3a$VCV[, 4] + 1 + pi^2/3)
mean(IC.b2);HPDinterval(IC.b2)


### RESCALE ESTIMATES BASED ON FIXED RESIDUAL VARIANCE ###

# rescale estimates by the estimated residual variance in order to obtain the
# posterior distributions of the parameters under the assumption that the 
# actual residual variance is equal to some other value.
c2 <- ((16 * sqrt(3))/(15 * pi))^2 # constant for logit link, c = 1 for probit

Sol.s <- m.3a$Sol/sqrt(1 + c2 * 1) # var = 1
Sol.s <- m.3a$Sol/sqrt(1 + c2 * 0) # var = 0

# Sol.s <- inv.logit(m.1$Sol/sqrt(1 + c2 * m.1$VCV[, 5]))
head(colnames(m.1$Sol))
mean(Sol.s[,"traitb1.1"]);HPDinterval(Sol.s[,"traitb1.1"])
mean(Sol.s[,"traitb2.1"]);HPDinterval(Sol.s[,"traitb2.1"])

### END ####

## NOTES ####

## UV PRIORS
# p <- list(G = list(G1=list(V=1, nu=0.002)), # Inv Wish for random intercepts (MCMCglmm default)
#           R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)), # parameter expanded G priors increase ESS of intercept and phylogenetic variance
#           R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1)), # reducing nu reduces ESS and gives inflated estimates of phylogenetic variance
#           R = list(V = 1, fix = 1))
# p <-  list(B = list(mu = c(0, 0), V = diag(2) * (1 + pi^2/3)), # prior to use when incomplete separation is a problem for fixed effects
#                    R = list(V = 1, fix = 1))
# p <- list(G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)), # parameter expanded prior recommended for binary data in MCMCMglmm course notes reduces ESS and gives inflated estimates of phylogenetic variance
#           R = list(V = 1,  fix = 1))


# p <- list(R = list(R1=list(V = diag(1), nu=0.002)),
#           G = list(G1=list(V = diag(1), nu=0.002)))
# m.1 <- MCMCglmm(g1 ~ g2, 
#                 random = ~animal, 
#                 rcov = ~units, 
#                 pedigree=t.toy,
#                 family = c("gaussian"), 
#                 nodes="ALL", data = d.toy, prior=p,
#                 nitt=35000, burnin=10000, thin=30, 
#                 pr=T, saveX = T, saveZ = T,
#                 verbose = F) 
# summary(m.1)
# 
# head(d.toy)
# # conditional mean
# W.1<-cBind(m.1$X, m.1$Z) # note X and Z are sparse so use cBind
# prediction.1<-W.1%*%posterior.mode(m.1$Sol)
# plot(d.toy$g1,prediction.1)
# plot(order(as.factor(d.toy$animal)), prediction.1@x)
# 
# # fixed effects only
# prediction.2<-m.1$X%*%posterior.mode(m.1$Sol[,1:2])
# plot(d.toy$g1,prediction.2);abline(0,1)
# plot(order(as.factor(d.toy$animal)), prediction.2@x)
# plot(order(as.factor(d.toy$animal)), d.toy$g2)
# 
# # cond mean has higher variance due to consideration of group level effects
# plot(prediction.2,prediction.1)
# var(prediction.1@x);var(prediction.2@x)
# 
# prediction.2@x
# 
# # group level effects
# plot(order(as.factor(d.toy$animal)),prediction.1-prediction.2);abline(0,0)
# 
# 
# ## using predict()
# #
# dummy <- data.frame(g2 = seq(0,1, length.out = 10),
#                     animal = paste0("s",rep(1:250, each = 10)),
#                     y=0)
# 
# # conditional on group level effects
# g1.c <- predict(m.1, dummy, marginal = NULL)
# # marginalizing over group level effects
# g1.m <- predict(m.1, dummy)
# 
# df.pred <- data.frame(animal=dummy$animal, g2=dummy$g2, g1.c, g1.m)
# 
# lattice::xyplot(g1.c + g1.m ~ g2 | animal, data = df.pred)[1:10]
# 
# head(df.pred)
# head(d.toy)
summary(m.1)

## END ####

#--------------------------------------------------------------#


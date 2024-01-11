library(tidyverse)
library(dplyr)
library(tidyr)
library(brms)
library(bayesplot)
library(ggpubr)
library(purrr)



rm(list = ls())
#fit <- readRDS("PMM_brms_2022_02_08.rds")

# choose model
fit <- fit_MR_2

# time to fit (seconds): slowest chain of four
rstan::get_elapsed_time(fit$fit) %>% apply(1,sum) %>% max

# names of model parameters, excluding random effects and log-probability
pars <- fit %>% as_tibble %>% select(-starts_with("r_"),-"lp__") %>% names
r_pars <- fit %>% as_tibble %>% select(starts_with("r_")) %>% names


#------------------------
# functions
#------------------------

# specify dplyr select function
select <- dplyr::select

# calculate coverage
get_coverage <- function(data) data %>%
  transmute(cov90 = (hh >= y_obs) & (ll <= y_obs),
            cov50 = (h >= y_obs) & (l <= y_obs)) %>% 
  map_dbl(sum) %>% `/`(nrow(data)) %>% `*`(100) %>% round(.,1)

# plot predictions
plot_pred <- function(x, y_obs = NULL, cov, ylab = NULL, xlab = NULL, title = NULL, subtitle = FALSE, arrange = FALSE) {
  if(!is.null(y_obs)) x <- x %>% mutate(y_obs = y_obs)
  if(arrange == TRUE) x <- x %>% arrange(m)
  x %>%  
  mutate(x = 1:n()) %>% 
  ggplot(aes(x)) +
  geom_linerange(aes(ymin = ll, ymax = hh), size = 1, col = blues9[3]) +
  geom_linerange(aes(ymin = l, ymax = h), size = 1, col = blues9[5]) +
    geom_point(aes(y=m), size = 2, col = blues9[7]) +
  geom_point(aes(y = y_obs), size = 2) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 10)) +
  labs(x = xlab,
       y = ylab,
       title = title,
       subtitle = ifelse(subtitle == TRUE, paste("Coverage: ", cov[1], "% at 90th quantile, ",cov[2],"% at 50th quantile"), ""))
}

#------------------------
# convergence diagnostics
#------------------------
mcmc_trace(fit, pars = pars)
mcmc_rhat(rhat(fit, pars = r_pars))
mcmc_neff(neff_ratio(fit, pars = r_pars))

rstan::check_hmc_diagnostics(fit$fit)

#----------------------------
# posterior predictive checks
#----------------------------

# CONDITIONAL - prediction conditional on species-level random-intercept estimates
ppC1 <- pp_check(fit, resp= "Nmass", type = "intervals") + labs(subtitle = "Response: N_mass")
ppC2 <- pp_check(fit, resp= "d13C", type = "intervals") + labs(subtitle = "Response: d13C")
ppC3 <- pp_check(fit, resp= "LMA", type = "intervals") + labs(subtitle = "Response: LMA")

covC1 <- get_coverage(ppC1$data)
covC2 <- get_coverage(ppC2$data)
covC3 <- get_coverage(ppC3$data)

ggarrange(plot_pred(ppC1$data %>% select(6:10), ppC1$data$y_obs, covC1, "N_mass", arrange = TRUE),
          plot_pred(ppC2$data %>% select(6:10), ppC2$data$y_obs, covC2, "d13C", arrange = TRUE),
          plot_pred(ppC3$data %>% select(6:10), ppC3$data$y_obs, covC3, "LMA", arrange = TRUE),
          ncol = 1)

# MARGINAL - prediction marginalised over Gaussian distribution of phylogenetic random effects
ppM1 <- pp_check(fit, resp= "Nmass", type = "intervals", allow_new_levels = T, sample_new_levels = "gaussian", newdata = fit$data %>% mutate(phylo = NA))
ppM2 <- pp_check(fit, resp= "d13C", type = "intervals", allow_new_levels = T, sample_new_levels = "gaussian", newdata = fit$data %>% mutate(phylo = NA))
ppM3 <- pp_check(fit, resp= "LMA", type = "intervals", allow_new_levels = T, sample_new_levels = "gaussian", newdata = fit$data %>% mutate(phylo = NA))

covM1 <- get_coverage(ppM1$data)
covM2 <- get_coverage(ppM2$data)
covM3 <- get_coverage(ppM3$data)

# arrange by predictive mean of the conditional model
ppM1$data <- ppM1$data %>% arrange(match(y_id, ppC1$data %>% arrange(m) %>% pull(y_id)))
ppM2$data <- ppM2$data %>% arrange(match(y_id, ppC2$data %>% arrange(m) %>% pull(y_id)))
ppM3$data <- ppM3$data %>% arrange(match(y_id, ppC3$data %>% arrange(m) %>% pull(y_id)))

ggarrange(plot_pred(ppM1$data %>% select(6:10), ppM1$data$y_obs, covM1, "N_mass"),
          plot_pred(ppM2$data %>% select(6:10), ppM2$data$y_obs, covM2, "d13C"),
          plot_pred(ppM3$data %>% select(6:10), ppM3$data$y_obs, covM3, "LMA"),
          ncol = 1)

# FIXED EFFECTS ONLY - (i.e., ignore random effects by setting them to zero)
ppF1 <- pp_check(fit, resp= "Nmass", type = "intervals", re_formula = NA) + labs(subtitle = "Response: N_mass")
ppF2 <- pp_check(fit, resp= "d13C", type = "intervals", re_formula = NA) + labs(subtitle = "Response: d13C")
ppF3 <- pp_check(fit, resp= "LMA", type = "intervals", re_formula = NA) + labs(subtitle = "Response: LMA")

covF1 <- get_coverage(ppF1$data)
covF2 <- get_coverage(ppF2$data)
covF3 <- get_coverage(ppF3$data)

# arrange by predictive mean of the conditional model
ppF1$data <- ppF1$data %>% arrange(match(y_id, ppC1$data %>% arrange(m) %>% pull(y_id)))
ppF2$data <- ppF2$data %>% arrange(match(y_id, ppC2$data %>% arrange(m) %>% pull(y_id)))
ppF3$data <- ppF3$data %>% arrange(match(y_id, ppC3$data %>% arrange(m) %>% pull(y_id)))

ggarrange(plot_pred(ppF1$data %>% select(6:10), ppF1$data$y_obs, covF1, "N_mass"),
          plot_pred(ppF2$data %>% select(6:10), ppF2$data$y_obs, covF2, "d13C"),
          plot_pred(ppF3$data %>% select(6:10), ppF3$data$y_obs, covF3, "LMA"),
          ncol = 1)


#----------------
# loo predictions
#----------------

# draw posterior samples
ppred <- posterior_predict(fit)
ppred1 <- ppred[,,1] # N_mass
ppred2 <- ppred[,,2] # d13C
ppred3 <- ppred[,,3] # LMA

# compute importance sampling weights
log_ratios1 <- -1*log_lik(fit, resp = "Nmass")
log_ratios2 <- -1*log_lik(fit, resp = "d13C")
log_ratios3 <- -1*log_lik(fit, resp = "LMA")
r_eff1 <- loo::relative_eff(exp(-log_ratios1), chain_id = rep(1:4, each = 1000))
r_eff2 <- loo::relative_eff(exp(-log_ratios2), chain_id = rep(1:4, each = 1000))
r_eff3 <- loo::relative_eff(exp(-log_ratios3), chain_id = rep(1:4, each = 1000))
psis1 <- loo::psis(log_ratios1, cores = 10, r_eff = r_eff1)
psis2 <- loo::psis(log_ratios2, cores = 10, r_eff = r_eff2)
psis3 <- loo::psis(log_ratios3, cores = 10, r_eff = r_eff3)

# generate loo quantile predictions via PSIS
loo_pred1 <- loo::E_loo(ppred1,psis_object = psis1, log_ratios = log_ratios1, type = "quantile", 
                        probs = c(0.05,0.25,0.5,0.75,0.95))
loo_pred2 <- loo::E_loo(ppred2,psis_object = psis2, log_ratios = log_ratios2, type = "quantile", 
                        probs = c(0.05,0.25,0.5,0.75,0.95))
loo_pred3 <- loo::E_loo(ppred3,psis_object = psis3, log_ratios = log_ratios3, type = "quantile", 
                        probs = c(0.05,0.25,0.5,0.75,0.95))

# collate loo estimates
loo_pred_probs1 <- loo_pred1$value %>% t %>% as.data.frame() %>% as_tibble()
loo_pred_probs2 <- loo_pred2$value %>% t %>% as.data.frame() %>% as_tibble()
loo_pred_probs3 <- loo_pred3$value %>% t %>% as.data.frame() %>% as_tibble()
colnames(loo_pred_probs1) <- colnames(loo_pred_probs2) <- colnames(loo_pred_probs3) <- c("ll","l","m","h","hh")

# compute coverage
cov1 <- loo_pred_probs1 %>% mutate(x = 1:n(), k = loo_pred1$pareto_k, y_obs = fit$data$N_mass) %>% 
  get_coverage() * 100/length(fit$data$animal); cov1 <- round(cov1)
cov2 <- loo_pred_probs2 %>% mutate(x = 1:n(), k = loo_pred2$pareto_k, y_obs = fit$data$d13C) %>% 
  get_coverage() * 100/length(fit$data$animal); cov2 <- round(cov2)
cov3 <- loo_pred_probs3 %>% mutate(x = 1:n(), k = loo_pred3$pareto_k, y_obs = fit$data$LMA) %>% 
  get_coverage() * 100/length(fit$data$animal); cov3 <- round(cov3)


# main plot
ggarrange(plot_pred(loo_pred_probs1, fit$data$N_mass, cov1, "N_mass", arrange = TRUE),
          plot_pred(loo_pred_probs2, fit$data$d13C, cov2, "d13C", arrange = TRUE),
          plot_pred(loo_pred_probs3, fit$data$LMA, cov3, "LMA", arrange = TRUE),
          ncol = 1)


#------------------
# residual analysis
#------------------
pe1 <- predictive_error(fit, resp = "Nmass") %>% as.data.frame() %>% map_dbl(mean)
pe2 <- predictive_error(fit, resp = "d13C") %>% as.data.frame() %>% map_dbl(mean)
pe3 <- predictive_error(fit, resp = "LMA") %>% as.data.frame() %>% map_dbl(mean)

tibble(r1 = pe1, r2 = pe2, r3 = pe3) %>% 
  ggplot() + 
  stat_qq(aes(sample = r1, col = "N_mass")) +
  stat_qq(aes(sample = r2, col = "d13C")) +
  stat_qq(aes(sample = r3, col = "LMA")) +
  scale_colour_manual(values = c(blues9[4], blues9[6], blues9[8])) +
  labs(col = "Response", subtitle = "QQ plots") +
  theme_classic()


# #------------------------------------
# # recovery of simulation parameters
# #------------------------------------
# x <- fit %>% as_tibble %>% select(-starts_with("r_"),-"lp__")
# true_values1 <- c(b_y1_Intercept = 1,
#                  b_y2_Intercept = 2,
#                  sd_animal__y1_Intercept = 1,
#                  sd_animal__y2_Intercept = 2, # sqrt needed here?
#                  cor_animal__y1_Intercept__y2_Intercept = 0.75,
#                  sigma_y1 = 2, # sqrt needed here?
#                  sigma_y2 = 1,
#                  rescor__y1__y2 = 0.25
#                  )
# 
# true_values2 <- c(b_y1_Intercept = 1,
#                  b_y2_Intercept = 2,
#                  sd_animal__y1_Intercept = 2,
#                  sd_animal__y2_Intercept = 4, # sqrt needed here?
#                  cor_animal__y1_Intercept__y2_Intercept = 0.75,
#                  sigma_y1 = 1, # sqrt needed here?
#                  sigma_y2 = 2,
#                  rescor__y1__y2 = 0.25
# )
# identical(names(true_values),names(x))
# mcmc_recover_hist(x, true_values2)
# 
# ## Condititonal prediction etc.
# fit
# data <- fit$data
# A <- fit$data2$A
# 
# # leave out first obs (manually)
# data_22 <- data[-1,]
# A_22 <- A[-1,-1]
# 
# fit_22 <- update(fit, newdata = data_22, data2 = list(A = A_22), cores = 4, recompile = T)
# saveRDS(fit_22, "fit_22.rds")
# fit_22
# 
# A_22_inv <- solve(A_22)
# 
# kronecker(A,diag(nrow = 2))
# 
# A_eig <- eigen(A)
# 
# Lambda <- A_eig$values %>% diag
# Lamba_inv <- (1/A_eig$values) %>% diag
# 
# Q <- A_eig$vectors
# 
# (Q) %*% Lambda %*% t(Q)  # eigendecomposition
# (Q) %*% Lambda_Inv %*% t(Q)
# 
# A <- matrix(c(1,3,3,5),2,2); A
# Sigma_phy <- matrix(c(1,0.5,0.5,1),2,2); Sigma_phy
# kronecker(t(A[1,]), Sigma_phy)
# kronecker(Sigma_phy,A)

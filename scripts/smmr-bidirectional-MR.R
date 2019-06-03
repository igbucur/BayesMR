##' Examine the robustness of the method with regard to IV validity

# 0. Load packages and functions ------------------------------------------

library(rstan)
library(bayesplot)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(configr)
source('utils/derive_sufficient_statistics.R')

gen_data_miv_sem <- function(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, seed = NULL) {
  
  set.seed(seed)
  
  G <- sapply(theta, function(t) rbinom(N, n, t)) # genetic variant
  U <- rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + rnorm(N, sd = sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + rnorm(N, sd = sigma_Y)
  Z <- cbind(1, G, X, Y)
  
  list(data = Z, SS = t(Z) %*% Z / N)
}

smmr_run_miv_sim <- function(SS, N, J, n, theta, mu_X = 0, mu_Y = 0, v_spike = 0.01, v_slab = 1, iter = 20000, chains = 10) {
  
  library(rstan)
  stanfit <- stan("src/multiple_IV.stan", data = list(
    N = N,
    J = J,
    n = n,
    v_spike = v_spike,
    v_slab = v_slab,
    SS = SS,
    mu_X = mu_X,
    mu_Y = mu_Y,
    theta = theta
  ), iter = iter, chains = chains, control = list(adapt_delta = 0.99, max_treedepth = 20))
  
  stanfit
}

write_polychord_config <- function(
  file_path, output_root,
  nlive, n_dims
) {
  
  config <- list()
  config[["algorithm settings"]] <- list(
    nlive = nlive,
    num_repeats = 2 * n_dims,
    do_clustering = "T",
    grade_frac = "1",
    precision_criterion = "0.001"
  )
  
  config[["posterior settings"]] <- list(
    posteriors = "T",
    equals = "T",
    cluster_posteriors = "F",
    boost_posterior = "5.0"
  )
  
  config[["output settings"]] <- list(
    base_dir = "/home/gbucur/dev/chains",
    file_root = output_root,
    write_resume = "F",
    write_dead = "F",
    feedback = 1,
    write_paramnames = "T"
  )
  
  config[["derived parameter settings"]] <- list()
  
  config[["prior settings"]] <- list()
  
  write.config(config, file_path)
  
  gen_prior_string <- function(idx) {
    paste0("P : p", idx, " | \\theta_{", idx, "} | 1 | uniform | 1 | 0.0 1.0")
  }
  
  for (i in 1:n_dims) {
    write(gen_prior_string(i), file_path, append = T)
  }
}

write_model_config <- function(
  file_path,
  J, SS_filename, sigmaG_filename,
  slab_precision = 1.0,
  spike_precision = 100.0,
  n_samples = 10000
) {
  
  config <- list()
  config[["model settings"]] <- list(
    SS_filename = SS_filename,
    sigmaG_filename = sigmaG_filename,
    slab_precision = slab_precision,
    spike_precision = spike_precision,
    instruments = J,
    observations = n_samples,
    model = 0
  )
  
  write.config(config, file_path)
}

biMR_save_data <- function(
  secondary_effect = 0, 
  file_root = sprintf('dev/PolyChord/data/bidirectional_alpha=%+1.2f', secondary_effect),
  seed = NULL) {
  
  J <- 2 # total number of generated IVs
  n <- 2 # binomial trials
  N <- 10000 # number of observations
  theta <- c(0.5, 0.5) # probability of exotic gene
  gamma <- c(1, secondary_effect) # instrument strengths
  alpha <- c(secondary_effect, 1) # pleiotropic effects
  beta <- 1 # parameter of interest
  sigma_X <- sigma_Y <- 1 # intrinsic variability
  kappa_X <- kappa_Y <- 1 # confounding coefficients
  mu_X <- mu_Y <- 0 # expected values
  
  sigma_G <- sqrt(n * theta * (1 - theta))
  write.table(sigma_G, paste0(file_root, '_SiG.txt'), row.names = FALSE, col.names = FALSE)
  
  data <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, seed)
  
  SSc <- data$SS
  write.table(SSc, paste0(file_root, '_SSc.txt'), row.names = FALSE, col.names = FALSE)
  SSr <- SSc[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]
  write.table(SSr, paste0(file_root, '_SSr.txt'), row.names = FALSE, col.names = FALSE)

  list(data_path = file_root, data_obj = data)
}

biMR_dir_fixed_params <- function(
  polychord_config_path, model_config_path, data_root, n_dims, n_live = 1000
) {

  write_polychord_config(polychord_config_path, paste0("biMR_PC_dir_sp2_", n_live, "_", basename(data_root)), n_live, n_dims)
  write_model_config(model_config_path, J, paste0(data_root, '_SSc.txt'), paste0(data_root, '_SiG.txt'))
  
  system(paste('dev/PolyChord/bin/polychord_MR_dir', polychord_config_path))
}

biMR_rev_fixed_params <- function(
  polychord_config_path, model_config_path, data_root, n_dims, n_live = 1000
) {
  
  write_polychord_config(polychord_config_path, paste0("biMR_PC_rev_sp2_", n_live, "_", basename(data_root)), n_live, n_dims)
  write_model_config(model_config_path, J, paste0(data_root, '_SSr.txt'), paste0(data_root, '_SiG.txt'))
  
  system(paste('dev/PolyChord/bin/polychord_MR_rev', polychord_config_path))
}

for (se in seq(-0.5, -0.1, 0.1)) {
  data_root <- biMR_save_data(se, seed = 1930)$data_path
  biMR_dir_fixed_params("dev/PolyChord/ini/config_dir.ini", "dev/PolyChord/ini/model_dir.ini", data_root, 11, 2000)
  biMR_rev_fixed_params("dev/PolyChord/ini/config_rev.ini", "dev/PolyChord/ini/model_rev.ini", data_root, 11, 2000)
}
# 


biMR_summary <- data.frame(
  secondary_effect= numeric(0),
  dir_gamma_hat = numeric(0), dir_gamma_se = numeric(0),
  dir_Gamma_hat = numeric(0), dir_Gamma_se = numeric(0),
  rev_gamma_hat = numeric(0), rev_gamma_se = numeric(0),
  rev_Gamma_hat = numeric(0), rev_Gamma_se = numeric(0),
  prob_dir = numeric(0)
)

biMR_data <- readRDS('results/bidirectional_MR_data.rds')
for (se in seq(-0.5, 0.5, 0.1)) {
  # biMR <- biMR_save_data(se, seed = 1930)
  biMR <- biMR_data[[as.character(se)]]

  lm_GX <- summary(lm(biMR$data_obj$data[, 4] ~ biMR$data_obj$data[, 2]))
  dir_gamma_hat <- lm_GX$coefficients[2, 1]
  dir_gamma_se <- lm_GX$coefficients[2, 2]

  lm_GY <- summary(lm(biMR$data_obj$data[, 5] ~ biMR$data_obj$data[, 2]))
  dir_Gamma_hat <- lm_GY$coefficients[2, 1]
  dir_Gamma_se <- lm_GY$coefficients[2, 2]

  lm_GX <- summary(lm(biMR$data_obj$data[, 4] ~ biMR$data_obj$data[, 3]))
  rev_Gamma_hat <- lm_GX$coefficients[2, 1]
  rev_Gamma_se <- lm_GX$coefficients[2, 2]

  lm_GY <- summary(lm(biMR$data_obj$data[, 5] ~ biMR$data_obj$data[, 3]))
  rev_gamma_hat <- lm_GY$coefficients[2, 1]
  rev_gamma_se <- lm_GY$coefficients[2, 2]
  
  Bf <- derive_polychord_Bayes_factor(
    sprintf('chains/bidirectional/direct/biMR_PC_dir_sp2_2000_bidirectional_alpha=%+1.2f.stats', se),
    sprintf('chains/bidirectional/reverse/biMR_PC_rev_sp2_2000_bidirectional_alpha=%+1.2f.stats', se)
  )

  biMR_summary[nrow(biMR_summary) + 1, ] <- c(
    secondary_effect = se,
    dir_gamma_hat = dir_gamma_hat,
    dir_gamma_se = dir_gamma_se,
    dir_Gamma_hat = dir_Gamma_hat,
    dir_Gamma_se = dir_Gamma_se,
    rev_gamma_hat = rev_gamma_hat,
    rev_gamma_se = rev_gamma_se,
    rev_Gamma_hat = rev_Gamma_hat,
    rev_Gamma_se = rev_Gamma_se,
    prob_dir = Bf$prob_dir
  )
}

library(dplyr)
biMR_summary <- biMR_summary %>%
  mutate(dir_beta_hat = dir_Gamma_hat / dir_gamma_hat, rev_beta_hat = rev_Gamma_hat / rev_gamma_hat) %>%
  mutate(dir_beta_se = dir_Gamma_se / dir_gamma_hat, rev_beta_se = rev_Gamma_se / rev_gamma_hat) %>%
  mutate(dir_beta_lower = dir_beta_hat - 2 * dir_beta_se, dir_beta_upper = dir_beta_hat + 2 * dir_beta_se) %>%
  mutate(rev_beta_lower = rev_beta_hat - 2 * rev_beta_se, rev_beta_upper = rev_beta_hat + 2 * rev_beta_se)

library(Hmisc)
latex(biMR_summary %>%
  filter(abs(secondary_effect) <= 0.5) %>%
  select(secondary_effect, dir_beta_hat, dir_beta_se, rev_beta_hat, rev_beta_se, prob_dir), rowname = NULL, file = "", digits = 3,
  colheads = c("Secondary effect", "$\\hat{\\ce}^{IV}_{X \\to Y}$", "$\\hat{\\sigma}^{IV}_{X \\to Y}$",
               "$\\hat{\\ce}^{IV}_{Y \\to X}$", "$\\hat{\\sigma}^{IV}_{Y \\to X}$", "$\\hat{p}(\\model_{X \\to Y} \\given \\data)$"))

# write_polychord_config("dev/PolyChord/ini/config_dir.ini", "MR", 1000, 9)
# write_model_config("dev/PolyChord/ini/model_dir.ini", 1, "expected_direction_SSc.txt", "expected_direction_SiG.txt")
# 
# sigma_G <- sqrt(n * theta * (1 - theta))
# write.table(sigma_G, 'dev/PolyChord/bidirectional_SiG.txt', row.names = FALSE, col.names = FALSE)

# TODO: see what happens with balanced pleiotropy
# see what happens when the correlation 

source('utils/PolyChord_interface.R')



       
# 1. Correct Direction ----------------------------------------------------

# parm <- get_matrix_SEM(J, n, theta, gamma, alpha, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
# SSc <- get_moments_matrix(J, parm$mu, parm$B, parm$C)

# SSc <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, 1930)
# write.table(SSc, 'dev/PolyChord/bidirectional_SSc.txt', row.names = FALSE, col.names = FALSE)
# 
# fitc <- smmr_run_miv_sim(SSc$SS, N, J, n, array(theta), chains = 2, iter = 10000)
# No pleiotropy (1, 0)
# -0.339955951051760E+005 +/-   0.199633068774385E+000
# -0.339959919749074E+005 +/-   0.153526208571663E+000
# Full pleiotropy (1, 1)
# PC: -0.339966851602529E+005 +/-   0.261851101453383E+000
# PC: -0.339967888696120E+005 +/-   0.151718039322212E+000

# levc <- c(-33910.365131782346, -33909.390909256501, -33909.051717157417, -33909.551850867661, -33911.082196318021)
# stdc <- c(0.26522881216627342, 0.26170274037905733, 0.26031196723177530, 0.26231136573092168, 0.26815730891358391) 
# 
# covc <- SSc[-1, -1] - SSc[-1, 1, drop = FALSE] %*% SSc[1, -1, drop = FALSE]

# 2. Reverse Direction ----------------------------------------------------
# SSr <- SSc[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]
# write.table(SSr, 'dev/PolyChord/bidirectional_SSr.txt', row.names = FALSE, col.names = FALSE)
# 
# fitr <- smmr_run_miv_sim(SSr, N, J, n, theta, chains = 2, iter = 5000)
# 
# levr <- c(-33911.872273660694, -33910.814989402716, -33909.804267330677, -33908.959235863455, -33909.283811595276)
# stdr <- c(0.27124940958613442, 0.26688342919799984, 0.26293213143543903, 0.26036808960990987, 0.26168343772207536)
# 
# covr <- SSr[-1, -1] - SSr[-1, 1, drop = FALSE] %*% SSr[1, -1, drop = FALSE]
# No pleiotropy (0, 1)
# -0.339975511105062E+005 +/-   0.222115699480298E+000
# Full pleiotropy (1, 1)
# PC: -0.339949890679235E+005 +/-   0.347869606391938E+000
# PC: -0.339951952968830E+005 +/-   0.151476980583303E+000

# df <- data.frame(cbind(levc, stdc, levr, stdr, a = c(0, 0.25, 0.5, 0.75, 1)))
# df <- df %>% mutate(lb = exp(levc - stdc - (levr + stdr)),
#               ub = exp(levc + stdc - (levr - stdr)),
#               bf = exp(levc - levr))
# 
# ggplot(df, aes(x=a, y= bf)) + 
#   geom_errorbar(aes(ymin = lb, ymax = ub), width=.1) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 1, col = "red") +
#   ylab("Correct versus reverse Bayes factor") +
#   xlab("Strength of pleiotropic effects") +
#   theme_tufte() +
#   theme(text = element_text(size = 15))
#   #geom_point()


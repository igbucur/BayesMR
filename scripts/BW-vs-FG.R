library(rstan)
source('utils/derive_sufficient_statistics.R')

# NOTE: We use only summary data
# MAF, IV-Exposure and IV-Outcome associations are from Greco et al. (2015).

# Exposure-Outcome association is from Daly et al. (2005).
# beta_hat <- -0.1555
beta_hat <- -0.00484
# beta_hat <- -0.06292 # -0.00484 - 3 * SD
# beta_hat <- -0.04356 # -0.00484 - 2 * SD
N <- 855

# Exposure-Outcome association is from Norris et al. (2011)
# beta_hat <- -0.025 
# N <- 6511

MAF <- c(0.61, 0.24, 0.5, 0.71, 0.7, 0.27, 0.74)
J <- 7
n <- 2

gamma_hat <- c(-0.072, -0.058, -0.045, -0.05, -0.024, -0.039, -0.037)
gamma_se <- c(0.007, 0.009, 0.007, 0.007, 0.009, 0.009, 0.010)
gamma_N <- c(34329, 34721, 41828, 42415, 26808, 29057, 23071)
Gamma_hat <- c(0.004, 0.024, 0.001, 0.010, 0.007, 0.007, 0.006)
Gamma_se <- c(0.004, 0.005, 0.004, 0.004, 0.004, 0.004, 0.005)
Gamma_N <- c(45727, 45726, 46186, 45056, 46171, 45062, 42074)

# Q Heterogeneity Test
# rma(Gamma_hat / gamma_hat, vi = (Gamma_hat / gamma_hat)^2 * ((Gamma_se / Gamma_hat)^2 + (gamma_se / gamma_hat)^2 - 2 * r * Gamma_se * gamma_se / (Gamma_hat * gamma_hat)), method = "FE")

SSS <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, N, MAF, n)
my_cov <- SSS[2:(J+3), 2:(J+3)] - SSS[2:(J+3), 1, drop = F] %*% SSS[1, 2:(J+3), drop = F]
my_cor <- cov2cor(my_cov)

# beta <- -0.1555
# # beta <- beta_hat
# var_G <- my_cov[1:J, 1:J]
# alpha <- Gamma_hat - beta * gamma_hat
# kappa_X <- kappa_Y <- c(sqrt(my_cov[J+1, J+2] - t(Gamma_hat) %*% var_G %*% gamma_hat - beta * my_cov[J+1, J+1]))
# sigma_X <- sqrt(my_cov[J+1, J+1] - t(gamma_hat) %*% var_G %*% gamma_hat - kappa_X^2)
# sigma_Y <- sqrt(my_cov[J+2, J+2] - t(Gamma_hat) %*% var_G %*% Gamma_hat - (beta * sigma_X)^2 - (kappa_Y + beta * kappa_X)^2)
# 
# set.seed(1915)
# nobs <- 1e6
# G <- sapply(MAF, rbinom, n = nobs, size = 2)
# U <- rnorm(nobs)
# 
# X <- G %*% gamma_hat + kappa_X * U + rnorm(nobs, sd = sigma_X)
# Y <- G %*% alpha + kappa_Y * U + beta * X + rnorm(nobs, sd = sigma_Y)
# Z <- cbind(1, G, X, Y)
# sum((t(Z) %*% Z / nobs - SSS)^2)

# Bridge sampling estimate of the log marginal likelihood: 1250.413 (sp4)
BW_FG_fit <- stan("src/multiple_IV.stan", data = list(
  N = N,
  J = J,
  n = n,
  v_spike = 1e-4,
  v_slab = 1,
  SS = SSS,
  mu_X = 0,
  mu_Y = 0,
  theta = MAF
), iter = 20000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))


# SSS2 <- MR_regression_coefficients_to_moments(1, gamma_hat[2], gamma_se[2], gamma_N[2], Gamma_hat[2], Gamma_se[2], Gamma_N[2], beta_hat, N, MAF[2], n)
# 
# BW_FG_fit_iv2 <- stan("src/multiple_IV.stan", data = list(
#   N = N,
#   J = 1,
#   n = n,
#   v_spike = 1e-2,
#   v_slab = 1,
#   SS = SSS2,
#   mu_X = 0,
#   mu_Y = 0,
#   theta = array(MAF[2])
# ), iter = 10000, chains = 1, control = list(adapt_delta = 0.99, max_treedepth = 20))

# reverse fit
# RSS <- SSS
# RSS[1:(J+1), c(J+2, J+3)] <- RSS[1:(J+1), c(J+3, J+2)]
# RSS[c(J+2, J+3), 1:(J+1)] <- RSS[c(J+3, J+2), 1:(J+1)]
# RSS[c(J+2, J+3), c(J+2, J+3)] <- RSS[c(J+3, J+2), c(J+3, J+2)]
# 
# BW_FG_fit_reverse <- stan("src/multiple_IV.stan", data = list(
#   N = N,
#   J = J,
#   n = n,
#   v_spike = 0.0001,
#   v_slab = 1,
#   SS = RSS,
#   mu_X = 0,
#   mu_Y = 0,
#   theta = MAF
# ), iter = 20000, chains = 5, control = list(adapt_delta = 0.99, max_treedepth = 20))

##' Examine the robustness of the method with regard to IV validity

# 0. Load packages and functions ------------------------------------------

library(rstan)
library(bayesplot)
source('utils/derive_sufficient_statistics.R')

gen_data_miv_sem <- function(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, seed = NULL) {
  
  set.seed(seed)
  
  G <- sapply(theta, function(t) rbinom(N, n, t)) # genetic variant
  U <- rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + rnorm(N, sd = sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + rnorm(N, sd = sigma_Y)
  Z <- cbind(1, G, X, Y)
  ESS <- t(Z) %*% Z / N # first and second-order moments
  
  ESS
}


gen_data_miv_sem_me <- function(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, tau_X = 1, tau_Y = 1, seed = NULL) {
  
  set.seed(seed)
  
  G <- sapply(theta, function(t) rbinom(N, n, t)) # genetic variant
  U <- rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + rnorm(N, sd = sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + rnorm(N, sd = sigma_Y)
  Xm <- X + rnorm(N, sd = tau_X)
  Ym <- Y + rnorm(N, sd = tau_Y)
  Zm <- cbind(1, G, Xm, Ym)
  ESSm <- t(Zm) %*% Zm / N # first and second-order moments
  
  ESSm
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


M <- 10 # total number of generated IVs
n <- 2 # binomial trials
N <- 100000 # number of observations

set.seed(1414)

theta <- runif(M, 0.0, 0.5) # probability of exotic gene
gamma <- rnorm(M) # instrument strengths
alpha <- rnorm(M, sd = 0.1) # pleiotropic effects
beta <- 1 # parameter of interest
kappa_X <- rnorm(1) # confounding coefficient on exposure
kappa_Y <- rnorm(1) # confounding coefficient on outcome
sigma_X <- abs(rnorm(1)) # intrinsic sd of exposure
sigma_Y <- abs(rnorm(1)) # intrinsic sd of outcome
# mu_X <- rnorm(1) # expected value of exposure
# mu_Y <- rnorm(1) # expected value of outcome

ESS <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)
ESSm <- gen_data_miv_sem_me(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)

fit <- smmr_run_miv_sim(ESS, N, M, n, array(theta), iter = 10000, chains = 3)
fitm <- smmr_run_miv_sim(ESSm, N, M, n, array(theta), iter = 10000, chains = 3)

# 1. Valid IV -------------------------------------------------------------

alpha_0 <- 0
param_0 <- get_matrix_SEM(M, n, theta, gamma, alpha_0, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
SS_0 <- get_moments_matrix(M, param_0$mu, param_0$B, param_0$C)

fit_0 <- smmr_run_miv_sim(SS_0, N, M, n, array(theta), iter = 10000, chains = 3, v_spike = 1e-4)

stanfit <- stan("src/posterior_direct.stan", data = list(
  N = N,
  n = n,
  v_spike = 0.0001,
  v_slab = 1,
  L = SS_0[1, 2],
  LL = SS_0[2, 2],
  LT = SS_0[2, 3:4],
  TT = SS_0[3:4, 3:4],
  theta = theta
), iter = 20000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))

fit_reverse <- stan("src/posterior_direct.stan", data = list(
  N = N,
  n = n,
  v_spike = 0.0001,
  v_slab = 1,
  L = SS_0[1, 2],
  LL = SS_0[2, 2],
  LT = SS_0[2, 4:3],
  TT = SS_0[4:3, 4:3],
  theta = theta
), iter = 20000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))



# 2. Weak Pleiotropy ------------------------------------------------------


# 3. Strong Pleiotropy ----------------------------------------------------



# 2. Simulation - 20/25 IVs are valid --------------------------------------
J <- 20 # number of valid IVs
alpha2 <- alpha; alpha2[seq(length.out = J)] <- 0

# ESS2 <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)
parm2 <- get_matrix_SEM(M, n, theta, gamma, alpha2, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
SS2 <- get_moments_matrix(M, parm2$mu, parm2$B, parm2$C)

fit2 <- smmr_run_miv_sim(SS2, N, M, n, theta)

# 3. Simulation - 15/25 IVs are valid --------------------------------------
J <- 15 # number of valid IVs
alpha3 <- alpha; alpha3[seq(length.out = J)] <- 0

# ESS2 <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)
parm3 <- get_matrix_SEM(M, n, theta, gamma, alpha3, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
SS3 <- get_moments_matrix(M, parm3$mu, parm3$B, parm3$C)

fit3 <- smmr_run_miv_sim(SS3, N, M, n, theta)


# 4. Simulation - 10/25 IVs are valid --------------------------------------
J <- 10 # number of valid IVs
alpha4 <- alpha; alpha4[seq(length.out = J)] <- 0

# ESS2 <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)
parm4 <- get_matrix_SEM(M, n, theta, gamma, alpha4, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
SS4 <- get_moments_matrix(M, parm4$mu, parm4$B, parm4$C)

fit4 <- smmr_run_miv_sim(SS4, N, M, n, theta)

# 5. Simulation - 05/25 IVs are valid --------------------------------------

J <- 5
alpha5 <- alpha; alpha5[seq(length.out = J)] <- 0

# ESS2 <- gen_data_miv_sem(N, n, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)
parm5 <- get_matrix_SEM(M, n, theta, gamma, alpha5, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
SS5 <- get_moments_matrix(M, parm5$mu, parm5$B, parm5$C)

fit5 <- smmr_run_miv_sim(SS5, N, M, n, theta)

# 6. Simulation - 00/25 IVs are valid --------------------------------------

parm <- get_matrix_SEM(M, n, theta, gamma, alpha, beta, 0, 0, kappa_X, kappa_Y, sigma_X, sigma_Y)
SS <- get_moments_matrix(M, parm$mu, parm$B, parm$C)

fit <- smmr_run_miv_sim(SS, N, M, n, theta)


# Rescale Parameters ------------------------------------------------------

sgamma <- gamma * sqrt(theta * (1 - theta) * n) / sigma_X
salpha <- alpha * sqrt(theta * (1 - theta) * n) / sigma_Y
sbeta <- beta * sigma_X / sigma_Y
skappa_X <- kappa_X / sigma_X
skappa_Y <- kappa_Y / sigma_Y


sigma <- matrix(c(sigma_X^2 + kappa_X^2, beta * (sigma_X^2 + kappa_X^2) + kappa_X * kappa_Y,
                  beta * (sigma_X^2 + kappa_X^2) + kappa_X * kappa_Y, sigma_Y^2 + (beta * sigma_X)^2 + (kappa_Y + beta * kappa_X)^2), 2, 2)


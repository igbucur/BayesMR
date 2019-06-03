#' @title Association between BMI and Parkinson's disease
#' @author Ioan Gabriel Bucur
#' @description Noyce, A. J. - Estimating the Causal Influence of Body Mass Index on Risk of Parkinson Disease

genetic_associations <- read.csv('data/SNP_bmi_parkinson.csv')


# choose 20 strongest instruments
# J <- 20

J <- nrow(genetic_associations)
dat <- genetic_associations[order(genetic_associations$gamma_hat, decreasing = T)[1:J],]

MAF <- dat$theta
gamma_hat <- dat$gamma_hat
gamma_se <- dat$gamma_se
Gamma_hat <- dat$Gamma_hat
Gamma_se <- dat$Gamma_se

beta_hat <- 0 # meta-analysis (Wang)
# beta_hat <- log(0.56) # Logroscino study
# beta_hat <- log(0.51) # Ma study

n <- 2
gamma_N <- 339224
Gamma_N <- 13708 + 95282
beta_N <- 430854 # meta-analysis
# beta_N <- 10812 # Logroscino study
# beta_N <- 425 # Ma study

N <- min(beta_N, gamma_N, Gamma_N)

# SSS <- matrix(0, J + 3, J + 3)
# 
# SSS[1, 1] <- 1
# SSS[1, 2:(J+1)] <- SSS[2:(J+1), 1] <- n * MAF
# SSS[2:(J+1), 2:(J+1)] <- n * n * MAF %*% t(MAF) + diag(n * MAF * (1 - MAF)) # E{GG^T}
# 
# SSS[1, J+2] <- SSS[J+2, 1] <- t(gamma_hat) %*% SSS[1, 2:(J+1)]
# SSS[1, J+3] <- SSS[J+3, 1] <- t(Gamma_hat) %*% SSS[1, 2:(J+1)]
# 
# SSS[J+2, 2:(J+1)] <- SSS[2:(J+1), J+2] <- SSS[2:(J+1), 2:(J+1)] %*% gamma_hat
# SSS[J+3, 2:(J+1)] <- SSS[2:(J+1), J+3] <- SSS[2:(J+1), 2:(J+1)] %*% Gamma_hat
# 
# SSS[J+2, J+2] <- mean((gamma_hat^2 + gamma_se^2 * gamma_N) * n * MAF * (1 - MAF)) + (n * t(gamma_hat) %*% MAF)^2
# SSS[J+2, J+3] <- SSS[J+3, J+2] <- mean((gamma_hat^2 + gamma_se^2 * gamma_N) * n * MAF * (1 - MAF)) * beta_hat + SSS[1, J+2] * SSS[J+3, 1]
# SSS[J+3, J+3] <- mean((Gamma_hat^2 + Gamma_se^2 * Gamma_N) * n * MAF * (1 - MAF)) + (n * t(Gamma_hat) %*% MAF)^2

SSS <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, beta_N, MAF, n)
#SSS[J+2, J+2] <- N / J * t(gamma_se) %*% diag(n * MAF * (1 - MAF)) %*% gamma_se + t(gamma_hat) %*% SSS[2:(J+1), 2:(J+1)] %*% gamma_hat
#SSS[J+2, J+3] <- SSS[J+3, J+2] <- N / J * t(gamma_se) %*% diag(n * MAF * (1 - MAF)) %*% gamma_se * beta_hat + beta_hat * t(gamma_hat) %*% diag(n * MAF * (1 - MAF)) %*% gamma_hat + SSS[1, J+2] * SSS[J+3, 1]
#SSS[J+3, J+3] <- N / J * t(Gamma_se) %*% diag(n * MAF * (1 - MAF)) %*% Gamma_se + t(Gamma_hat) %*% SSS[2:(J+1), 2:(J+1)] %*% Gamma_hat

library(rstan)
Pks_BMI_fit <- stan("src/multiple_IV.stan", data = list(
  N = N,
  J = J,
  n = n,
  v_spike = 1,
  v_slab = 1,
  SS = SSS,
  mu_X = 0,
  mu_Y = 0,
  theta = MAF
), iter = 20000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))


# library(parallel)
# library(rstan)
# 
# if (.Platform$OS.type == "unix") {
#   sflist <- 
#     mclapply(1:20, mc.cores = 20, 
#              function(i) {
#                stan("src/multiple_IV.stan", data = list(
#                  N = N,
#                  J = J,
#                  n = n,
#                  v_spike = 0.0001,
#                  v_slab = 1,
#                  SS = SSS,
#                  mu_X = 0,
#                  mu_Y = 0,
#                  theta = MAF
#                ), iter = 10000, chains = 1, chain_id = i, refresh = -1,
#                control = list(adapt_delta = 0.99, max_treedepth = 20))
#              }
#     )
#   fit_fake <- sflist2stanfit(sflist)
# }

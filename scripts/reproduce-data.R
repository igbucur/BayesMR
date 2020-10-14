
# 0. Setup ----------------------------------------------------------------



source('R/utils.R')

# 5.4 - Sensitivity to Nonlinearity and Nonnormality - Tables 2 and 3 -----


## Robustness to nonlinearity

# In the limit A -> infinity, the function becomes f(x) := x
parametric_nonlinearity <- function(x, A) {
  ifelse(A == Inf, x, A * tanh(x / A))
}

for (A in c(2^seq(-1, 5), Inf)) {
  set.seed(1918)
  
  samples <- 1e4
  MAF <- 0.3
  n <- 2
  J <- 1 # number of instruments
  G <- rbinom(samples, n, MAF)
  U <- rnorm(samples)
  X <- G + 0.1 * U + rnorm(samples)
  # we replace the linear Y-X relationship with a parametric nonlinearity
  Y <- 0.1 * G + parametric_nonlinearity(X, A) + 0.1 * U + rnorm(samples)
  
  Z <- cbind(1, G, X, Y)
  SSS <- t(Z) %*% Z / samples
  RSS <- SSS[c(1, 2, 4, 3), c(1, 2, 4, 3)]
  sigma_G <- sd(G) # sigma_G <- sqrt(MAF * (1 - MAF) * n)
  
  file_SS_expected <- sprintf('data/robustness_nonlinearity_A=%3.1f_SS_expected.txt', A)
  file_SS_reverse <- sprintf('data/robustness_nonlinearity_A=%3.1f_SS_reverse.txt', A)
  file_sigma_G <- sprintf('data/robustness_nonlinearity_A=%3.1f_sigma_G.txt', A)
  
  write.table(SSS, file = file_SS_expected, row.names = FALSE, col.names = FALSE)
  write.table(RSS, file = file_SS_reverse, row.names = FALSE, col.names = FALSE)
  write.table(sigma_G, file = file_sigma_G, row.names = FALSE, col.names = FALSE)
  
  # Write configuration file for expected direction model
  ini_expected <- sprintf('ini/robustness_nonlinearity_A=%3.1f_expected.ini', A)
  
  write_BayesMR_configuration_file(
    config_filename = ini_expected,
    SS_filename = file_SS_expected,
    sigmaG_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 500,
      num_repeats = 50
    )
  )
  # Sample from expected direction model using PolyChord
  run_system_command(paste('./BayesMR', ini_expected))
  
  
  # Write configuration file for reverse direction model
  ini_reverse <- sprintf('ini/robustness_nonlinearity_A=%3.1f_reverse.ini', A)
  
  write_BayesMR_configuration_file(
    config_filename = ini_reverse,
    SS_filename = file_SS_reverse,
    sigmaG_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 500,
      num_repeats = 50
    )
  )
  # Sample from reverse direction model using PolyChord
  run_system_command(paste('./BayesMR', ini_reverse))
}


## Robustness to Nonnormality

# NOTE: nu == Inf corresponds to the normal case

for (nu in c(2^seq(0, 6), Inf)) {
  
  set.seed(1918)
  
  samples <- 1e4
  MAF <- 0.3
  n <- 2
  J <- 1 # number of instruments
  G <- rbinom(samples, n, MAF)
  U <- rnorm(samples)
  X <- G + 0.1 * U + rnorm(samples)
  # We replace the normal noise for Y with t-distributed noise
  Y <- 0.1 * G + X + 0.1 * U + rt(samples, df = nu)
  
  Z <- cbind(1, G, X, Y)
  SSS <- t(Z) %*% Z / samples
  RSS <- SSS[c(1, 2, 4, 3), c(1, 2, 4, 3)]
  sigma_G <- sd(G) # sigma_G <- sqrt(MAF * (1 - MAF) * n)
  
  # sigma_G stays the same
  file_SS_expected <- sprintf('data/robustness_nonnormality_nu=%3.1f_SS_expected.txt', nu)
  file_SS_reverse <- sprintf('data/robustness_nonnormality_nu=%3.1f_SS_reverse.txt', nu)
  file_sigma_G <- sprintf('data/robustness_nonnormality_nu=%3.1f_sigma_G.txt', nu)
  
  write.table(SSS, file = file_SS_expected, row.names = FALSE, col.names = FALSE)
  write.table(RSS, file = file_SS_reverse, row.names = FALSE, col.names = FALSE)
  write.table(sigma_G, file = file_sigma_G, row.names = FALSE, col.names = FALSE)
  
  
  
  # Write configuration file for expected direction model
  ini_expected <- sprintf('ini/robustness_nonnormality_nu=%3.1f_expected.ini', nu)
  
  write_BayesMR_configuration_file(
    config_filename = ini_expected,
    SS_filename = file_SS_expected,
    sigmaG_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 500,
      num_repeats = 50
    )
  )
  # Sample from expected direction model using PolyChord
  run_system_command(paste('./BayesMR', ini_expected))
  
  
  # Write configuration file for reverse direction model
  ini_reverse <- sprintf('ini/robustness_nonnormality_nu=%3.1f_reverse.ini', nu)
  
  write_BayesMR_configuration_file(
    config_filename = ini_reverse,
    SS_filename = file_SS_reverse,
    sigmaG_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 500,
      num_repeats = 50
    )
  )
  # Sample from reverse direction model using PolyChord
  run_system_command(paste('./BayesMR', ini_reverse))
  
}

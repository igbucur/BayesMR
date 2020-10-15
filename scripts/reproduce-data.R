
# 0. Setup ----------------------------------------------------------------

source('R/utils.R')
source('R/PolyChord_interface.R')
source('R/write_BayesMR_configuration_file.R')

# WARNING: From a certain point on, there appeared issues with the PolyChordLite
# implementation whereby the samples were not saved properly at the end of the run,
# but remained in the temp files. That is why, we read samples from files
# ending with _temp_equal_weights.txt instead of _equal_weights.txt. This may
# not apply to all systems, so if [root]_temp_equal_weights.txt cannot be found,
# try changing it to [root]_equal_weights.txt.

if (!dir.exists("chains")) dir.create("chains/clusters", recursive = TRUE)

# 4.4 - Example: Near-conditional independence (near-LCD) - Figure 14 -----

set.seed(1918)

samples <- 1e4
MAF <- 0.3
n <- 2
J <- 1 # number of instruments
G <- rbinom(samples, n, MAF)
U <- rnorm(samples)
X <- G + 0.1 * U + rnorm(samples)
Y <- 0.1 * G + X + 0.1 * U + rnorm(samples)

Z <- cbind(1, G, X, Y)
SSS <- t(Z) %*% Z / samples
RSS <- SSS[c(1, 2, 4, 3), c(1, 2, 4, 3)]
sigma_G <- sd(G) # sigma_G <- sqrt(MAF * (1 - MAF) * n)

fileroot <- 'near_LCD_example'

file_SS_expected <- sprintf('data/%s_SS_expected.txt', fileroot)
file_SS_reverse <- sprintf('data/%s_SS_reverse.txt', fileroot)
file_sigma_G <- sprintf('data/%s_sigma_G.txt', fileroot)

write.table(SSS, file = file_SS_expected, row.names = FALSE, col.names = FALSE)
write.table(RSS, file = file_SS_reverse, row.names = FALSE, col.names = FALSE)
write.table(sigma_G, file = file_sigma_G, row.names = FALSE, col.names = FALSE)

ini_expected <- sprintf('ini/%s_expected.ini', fileroot)

write_BayesMR_configuration_file(
  config_filename = ini_expected,
  SS_filename = file_SS_expected,
  sigma_G_filename = file_sigma_G,
  num_instruments = 1,
  num_observations = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

ini_reverse <- sprintf('ini/%s_reverse.ini', fileroot)
write_BayesMR_configuration_file(
  config_filename = ini_reverse,
  SS_filename = file_SS_reverse,
  sigma_G_filename = file_sigma_G,
  num_instruments = 1,
  num_observations = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

run_system_command(paste("./BayesMR", ini_expected))
run_system_command(paste("./BayesMR", ini_reverse))

near_LCD_example <- list(
  samples_expected = read_PolyChord_samples(
    sprintf('chains/%s_expected_temp_equal_weights.txt', fileroot), quantile_transform = FALSE
  ),
  samples_reverse = read_PolyChord_samples(
    sprintf('chains/%s_reverse_temp_equal_weights.txt', fileroot), quantile_transform = FALSE
  ),
  evidence_comparison = derive_polychord_Bayes_factor(
    sprintf('chains/%s_expected.stats', fileroot), 
    sprintf('chains/%s_reverse.stats', fileroot)
  )
)

names(near_LCD_example$evidence_comparison) <- c(
  "Bayes_factor_expected_versus_reverse",
  "Bayes_factor_lower_standard_deviation",
  "Bayes_factor_upper_standard_deviation",
  "probability_expected",
  "probability_reverse"
)

save(near_LCD_example, file = 'data/near_LCD_example_reproduced.RData')

# 5.1 - Estimation Robustness - Figure 15 ---------------------------------

J <- 5 # total number of generated IVs
n <- 2 # binomial trials
N <- 100000 # number of observations

set.seed(1414)

theta <- runif(J, 0.0, 0.5) # probability of exotic gene
gamma <- rnorm(J, sd = 0.1) # instrument strengths
alpha <- rnorm(J, sd = 0.1) # pleiotropic effects
beta <- 1 # causal effect (parameter of interest)
kappa_X <- rnorm(1) # confounding coefficient on exposure
kappa_Y <- rnorm(1) # confounding coefficient on outcome
sigma_X <- abs(rnorm(1)) # intrinsic sd of exposure
sigma_Y <- abs(rnorm(1)) # intrinsic sd of outcome

estimation_robustness <- list()

for (K in seq(J, 0, -(J/5))) {
  alpha_K <- alpha
  alpha_K[1:K] <- 0
  
  fileroot <- sprintf('estimation_robustness_%d%%_valid', K / J * 100)
  
  # Sufficient statistics
  dat <- generate_data_BayesMR_model(N, n, theta, gamma, alpha_K, beta, kappa_X, kappa_Y, sigma_X, sigma_Y)
  SS_filename <- paste0("data/", fileroot, "_SS_expected.txt")
  write.table(dat$SS, SS_filename, row.names = FALSE, col.names = FALSE)
  
  # Estimated standard deviation of (independent) genetic variables G
  sigma_G <- sqrt((1 - dat$SS[1, 2:(J+1)] / n) * dat$SS[1, 2:(J+1)])
  sigma_G_filename <- paste0("data/", fileroot, "_sigma_G.txt")
  write.table(sigma_G, sigma_G_filename, row.names = FALSE, col.names = FALSE)
  
  ini_filename <- paste0("ini/", fileroot, ".ini")
  
  write_BayesMR_configuration_file(
    config_filename = ini_filename,
    SS_filename = SS_filename,
    sigma_G_filename = sigma_G_filename,
    num_instruments = J, num_observations = N
  )
  
  run_system_command(paste("./BayesMR", ini_filename))
  
  PolyChord_samples <- read_PolyChord_samples(
    sprintf('chains/%s_temp_equal_weights.txt', fileroot), quantile_transform = FALSE
  )
  estimation_robustness[[sprintf("%d%%", K / J * 100)]] <- 
    extract_beta_PolyChord_samples(PolyChord_samples)
}

save(estimation_robustness, file = 'data/estimation_robustness_reproduced.RData')

# 5.2 - Sensitivity to Prior Hyperparameters - Figure 16 ------------------

# Generate data for weak pleiotropy, weak confounding scenario
set.seed(1918)
samples <- 1e4
MAF <- 0.3
n <- 2
G <- rbinom(samples, n, MAF)
U <- rnorm(samples)
X <- G + 0.1 * U + rnorm(samples)
Y <- 0.1 * G + X + 0.1 * U + rnorm(samples)

Z <- cbind(1, G, X, Y)
sigma_G <- sqrt(MAF * (1 - MAF) * n)
SS_expected <- t(Z) %*% Z / samples
SS_reverse <- SS_expected[c(1, 2, 4, 3), c(1, 2, 4, 3)]

fileroot <- "sensitivity_to_hyperparameters"

SS_filename <- sprintf('data/%s_SS.txt', fileroot)
sigma_G_filename <- sprintf('data/%s_sigma_G.txt', fileroot)

write.table(SS_expected, file = SS_filename, row.names = FALSE, col.names = FALSE)
write.table(sigma_G, file = sigma_G_filename, row.names = FALSE, col.names = FALSE)

PolyChord_beta_samples <- list()

for (i in 0:6) {
  spike_precision <- 10^i
  
  ini_filename <- sprintf('ini/%s.ini', fileroot)
  write_BayesMR_configuration_file(
    config_filename = ini_filename, 
    SS_filename = SS_filename,
    sigma_G_filename = sigma_G_filename,
    num_instruments = 1,
    num_observations = N,
    spike_precision = spike_precision
  )
  
  run_system_command(paste('./BayesMR', ini_filename))
  
  PolyChord_beta_samples[[as.character(i)]] <- quick_derive_beta(
    utils::read.table(sprintf('chains/%s_temp_equal_weights.txt', fileroot))
  )
  
}

sensitivity_to_hyperparameters <- data.frame(spike = integer(0), samples = numeric(0))
for (i in 0:6) {
  sensitivity_to_hyperparameters <- rbind(
    sensitivity_to_hyperparameters, 
    cbind(rep(i, length(PolyChord_beta_samples[[as.character(i)]])), 
    PolyChord_beta_samples[[as.character(i)]])
  )
}
names(sensitivity_to_hyperparameters) <- c("spike", "beta")
sensitivity_to_hyperparameters$spike <- as.factor(sensitivity_to_hyperparameters$spike)

save(sensitivity_to_hyperparameters, file = 'data/sensitivity_to_hyperparameters_reproduced.RData')

# 5.3 - Model Averaging Robustness - Table 1 ------------------------------

model_averaging_robustness <- sapply(seq(-0.5, 0.5, 0.1), function(delta) {
  
  fileroot <- sprintf("model_averaging_robustness_delta=%+1.2f", delta)
  
  J <- 2 # total number of generated IVs
  n <- 2 # binomial trials
  N <- 10000 # number of observations
  theta <- c(0.5, 0.5) # probability of exotic gene
  gamma <- c(1, delta) # instrument strengths
  alpha <- c(delta, 1) # pleiotropic effects
  beta <- 1 # parameter of interest
  sigma_X <- sigma_Y <- 1 # intrinsic variability
  kappa_X <- kappa_Y <- 1 # confounding coefficients
  mu_X <- mu_Y <- 0 # expected values
  
  dat <- generate_data_BayesMR_model(N, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, 1930)
  SS_expected <- dat$SS
  SS_reverse <- dat$SS[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]

  SS_expected_filename <- sprintf('data/%s_SS_expected.txt', fileroot)
  SS_reverse_filename <- sprintf('data/%s_SS_reverse.txt', fileroot)
  sigma_G_filename <- sprintf('data/%s_sigma_G.txt', fileroot)
  
  write.table(sqrt(n * theta * (1 - theta)), sigma_G_filename, row.names = FALSE, col.names = FALSE)
  write.table(SS_expected, SS_expected_filename, row.names = FALSE, col.names = FALSE)
  write.table(SS_reverse, SS_reverse_filename, row.names = FALSE, col.names = FALSE)
  
  ini_expected_filename <- sprintf('ini/%s_expected.ini', fileroot)
  write_BayesMR_configuration_file(
    config_filename = ini_expected_filename, 
    SS_filename = SS_expected_filename,
    sigma_G_filename = sigma_G_filename,
    num_instruments = 2,
    num_observations = N
  )
  run_system_command(paste('./BayesMR', ini_expected_filename))
  
  ini_reverse_filename <- sprintf('ini/%s_reverse.ini', fileroot)
  write_BayesMR_configuration_file(
    config_filename = ini_reverse_filename, 
    SS_filename = SS_reverse_filename,
    sigma_G_filename = sigma_G_filename,
    num_instruments = 2,
    num_observations = N
  )
  run_system_command(paste('./BayesMR', ini_reverse_filename))
  
  lm_GX_expected <- summary(lm(dat$data[, 4] ~ dat$data[, 2]))
  expected_gamma_hat <- lm_GX_expected$coefficients[2, 1]
  expected_gamma_std_err <- lm_GX_expected$coefficients[2, 2]
  
  lm_GY_expected <- summary(lm(dat$data[, 5] ~ dat$data[, 2]))
  expected_Gamma_hat <- lm_GY_expected$coefficients[2, 1]
  expected_Gamma_std_err <- lm_GY_expected$coefficients[2, 2]
  
  lm_GX_reverse <- summary(lm(dat$data[, 4] ~ dat$data[, 3]))
  reverse_gamma_hat <- lm_GX_reverse$coefficients[2, 1]
  reverse_gamma_std_err <- lm_GX_reverse$coefficients[2, 2]
  
  lm_GY_reverse <- summary(lm(dat$data[, 5] ~ dat$data[, 3]))
  reverse_Gamma_hat <- lm_GY_reverse$coefficients[2, 1]
  reverse_Gamma_std_err <- lm_GY_reverse$coefficients[2, 2]
  
  Bayes_factor <- derive_polychord_Bayes_factor(
    sprintf('chains/%s_expected.stats', fileroot),
    sprintf('chains/%s_reverse.stats', fileroot)
  )
  
  c(
    delta = delta,
    dir_gamma_hat = dir_gamma_hat,
    dir_gamma_se = dir_gamma_se,
    dir_Gamma_hat = dir_Gamma_hat,
    dir_Gamma_se = dir_Gamma_se,
    rev_gamma_hat = rev_gamma_hat,
    rev_gamma_se = rev_gamma_se,
    rev_Gamma_hat = rev_Gamma_hat,
    rev_Gamma_se = rev_Gamma_se,
    probability_expected = Bayes_factor$prob_dir
  )
})

save(model_averaging_robustness, file = 'data/model_averaging_robustness_reproduced.RData')


# 5.4 - Sensitivity to Nonlinearity and Nonnormality - Tables 2 and 3 -----

## Robustness to nonlinearity
sensitivity_to_nonlinearity <- list()

# In the limit A -> infinity, the function becomes f(x) := x
parametric_nonlinearity <- function(x, A) {
  ifelse(A == Inf, x, A * tanh(x / A))
}

for (A in c(2^seq(-1, 5), Inf)) {
  
  # Simulate a weak pleiotropy, weak confounding scenario
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
    sigma_G_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 100,
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
    sigma_G_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 500,
      num_repeats = 50
    )
  )
  # Sample from reverse direction model using PolyChord
  run_system_command(paste('./BayesMR', ini_reverse))
  
  # Read and collect results from PolyChord output
  results <- list()
  results$expected$PolyChord_samples <- 
    read.table(sprintf('chains/robustness_nonlinearity_A=%3.1f_expected_temp_equal_weights.txt', A))
  results$expected$beta_posterior <- quick_derive_beta(results$expected$PolyChord_samples)
  results$reverse$PolyChord_samples <- 
    read.table(sprintf('chains/robustness_nonlinearity_A=%3.1f_reverse_temp_equal_weights.txt', A))
  results$direction_evidence <- derive_polychord_Bayes_factor(
    sprintf('chains/robustness_nonlinearity_A=%3.1f_expected.stats', A),
    sprintf('chains/robustness_nonlinearity_A=%3.1f_reverse.stats', A)
  )
  names(results$direction_evidence) <- c(
    'Bayes_factor', 'Bayes_factor_lower_error_bar', "Bayes_factor_upper_error_bar", 
    "probability_expected", "probability_reverse")
  
  results$direction_evidence$probability_expected_lower_error_bar <- 
    results$direction_evidence$Bayes_factor_lower_error_bar / 
    (1 + results$direction_evidence$Bayes_factor_lower_error_bar)
  
  results$direction_evidence$probability_expected_upper_error_bar <- 
    results$direction_evidence$Bayes_factor_upper_error_bar / 
    (1 + results$direction_evidence$Bayes_factor_upper_error_bar)
  
  sensitivity_to_nonlinearity[[paste0("A=", A)]] <- results
}


## Sensitivity to Nonnormality
sensitivity_to_nonnormality <- list()

# NOTE: nu == Inf corresponds to the normal case

for (nu in c(2^seq(0, 6), Inf)) {
  
  # Simulate a weak pleiotropy, weak confounding scenario
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
    sigma_G_filename = file_sigma_G,
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
    sigma_G_filename = file_sigma_G,
    num_instruments = J,
    num_observations = samples,
    PolyChord_control = list(
      num_live_points = 500,
      num_repeats = 50
    )
  )
  # Sample from reverse direction model using PolyChord
  run_system_command(paste('./BayesMR', ini_reverse))
  
  
  
  # Read and collect results from PolyChord output
  results <- list()
  results$expected$PolyChord_samples <- 
    read.table(sprintf('chains/robustness_nonnormality_nu=%3.1f_expected_temp_equal_weights.txt', nu))
  results$expected$beta_posterior <- quick_derive_beta(results$expected$PolyChord_samples)
  results$reverse$PolyChord_samples <- 
    read.table(sprintf('chains/robustness_nonnormality_nu=%3.1f_reverse_temp_equal_weights.txt', nu))
  results$direction_evidence <- derive_polychord_Bayes_factor(
    sprintf('chains/robustness_nonnormality_nu=%3.1f_expected.stats', nu),
    sprintf('chains/robustness_nonnormality_nu=%3.1f_reverse.stats', nu)
  )
  names(results$direction_evidence) <- c(
    'Bayes_factor', 'Bayes_factor_lower_error_bar', "Bayes_factor_upper_error_bar", 
    "probability_expected", "probability_reverse")
  
  results$direction_evidence$probability_expected_lower_error_bar <- 
    results$direction_evidence$Bayes_factor_lower_error_bar / 
    (1 + results$direction_evidence$Bayes_factor_lower_error_bar)
  
  results$direction_evidence$probability_expected_upper_error_bar <- 
    results$direction_evidence$Bayes_factor_upper_error_bar / 
    (1 + results$direction_evidence$Bayes_factor_upper_error_bar)
  
  sensitivity_to_nonnormality[[paste0("nu=", nu)]] <- results
}

save(sensitivity_to_nonlinearity, sensitivity_to_nonnormality, 
     file = 'data/sensitivity_to_parametric_assumptions_reproduced.RData')


# 6.2 - Effect of birth weight on fasting glucose - Figures 17 and 18 -----

# Extract BayesMR input from summary statistics

# Exposure-Outcome association is from Daly et al. (2005).
beta_hat <- -0.1555
N <- 855 # number of observations

EAF <- c(0.61, 0.24, 0.5, 0.71, 0.7, 0.27, 0.74) # effect allele frequencies
J <- 7 # number of instruments
n <- 2 # number of alleles (trials)

# EAF, IV-Exposure and IV-Outcome associations are from Greco et al. (2015).
gamma_hat <- c(-0.072, -0.058, -0.045, -0.05, -0.024, -0.039, -0.037) # G-X association
gamma_se <- c(0.007, 0.009, 0.007, 0.007, 0.009, 0.009, 0.010) # G-X std. err.
gamma_N <- c(34329, 34721, 41828, 42415, 26808, 29057, 23071) # G-X num. obs.
Gamma_hat <- c(0.004, 0.024, 0.001, 0.010, 0.007, 0.007, 0.006) # G-Y association
Gamma_se <- c(0.004, 0.005, 0.004, 0.004, 0.004, 0.004, 0.005) # G-Y std. err.
Gamma_N <- c(45727, 45726, 46186, 45056, 46171, 45062, 42074) # G-Y num. obs.

birth_weight_fasting_glucose_SS <- 
  MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, N, MAF, n)

# Write BayesMR input files

fileroot <- 'birth_weight_fasting_glucose'

ini_filename <- sprintf('ini/%s.ini', fileroot)
SS_filename <- sprintf('data/%s_SS.txt', fileroot)
sigma_G_filename <- sprintf('data/%s_sigma_G.txt', fileroot)

# Save first and second order statistics for the example
write.table(
  MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, N, MAF, n),
  file = SS_filename, row.names = FALSE, col.names = FALSE
)
write.table(
  sqrt(n * EAF * (1-EAF)), file = sigma_G_filename, row.names = FALSE, col.names = FALSE
)

# Posterior obtained with spike-and-slab prior (Figure 17) ----------------------------

write_BayesMR_configuration_file(
  config_filename = ini_filename,
  SS_filename = SS_filename ,
  sigma_G_filename = sigma_G_filename,
  num_instruments = 7,
  num_observations = 855,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 500)
)

run_system_command(paste('./BayesMR', ini_filename))

birth_weight_fasting_glucose_default <- 
  read_PolyChord_samples(sprintf('chains/%s_temp_equal_weights.txt', fileroot))

# Posterior obtained under IV assumptions (Figure 18) ---------------------------------

write_BayesMR_configuration_file(
  config_file_name = ini_filename,
  SS_filename = SS_expected_filename ,
  sigmaG_filename = sigma_G_filename,
  num_instruments = 7,
  num_observations = 855,
  spike_precision = 10000,
  model_type = 1, # use IV assumptions
  PolyChord_control = list(num_live_points = 500)
)

run_system_command(paste('./BayesMR', ini_filename))

birth_weight_fasting_glucose_IV <- 
  read_PolyChord_samples(sprintf('chains/%s_temp_equal_weights.txt', fileroot))

# Save data files
save(birth_weight_fasting_glucose_default, birth_weight_fasting_glucose_IV,
     file = 'data/birth_weight_fasting_glucose_reproduced.RData')


# 6.3 Effect of BMI on the risk of PD - Figures 19-21 ---------------------

# Extract BayesMR input from summary statistics

BMI_genetic_associations <- read.csv('inst/extdata/SNP_bmi_parkinson.csv')
J <- nrow(BMI_genetic_associations)
BMI_genetic_associations <- BMI_genetic_associations[order(BMI_genetic_associations$gamma_hat, decreasing = T)[1:J],]

EAF <- BMI_genetic_associations$theta
gamma_hat <- BMI_genetic_associations$gamma_hat
gamma_se <- BMI_genetic_associations$gamma_se
Gamma_hat <- BMI_genetic_associations$Gamma_hat
Gamma_se <- BMI_genetic_associations$Gamma_se
beta_hat <- 0 # meta-analysis (Wang et al. 2015)

n <- 2
gamma_N <- 339224
Gamma_N <- 13708 + 95282
beta_N <- 430854 # meta-analysis (Wang et al. 2015)

N <- min(beta_N, gamma_N, Gamma_N)

# Write BayesMR input files

fileroot <- 'Parkinson_BMI'

ini_filename <- sprintf('ini/%s.ini', fileroot)
SS_expected_filename <- sprintf('data/%s_SS_expected.txt', fileroot)
sigma_G_filename <- sprintf('data/%s_sigma_G.txt', fileroot)

SS_expected <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, beta_N, EAF, n)
write.table(SS_expected, file = SS_expected_filename, row.names = FALSE, col.names = FALSE)
write.table(sqrt(n * EAF * (1 - EAF)), file = sigma_G_filename, row.names = FALSE, col.names = FALSE)

# Writing the PolyChord configuration file

fileroot <- 'Parkinson_BMI'

ini_filename <- sprintf('ini/%s.ini', fileroot)
SS_expected_filename <- sprintf('data/%s_SS_expected.txt', fileroot)
sigma_G_filename <- sprintf('data/%s_sigma_G.txt', fileroot)


write_BayesMR_configuration_file(
  config_filename = ini_filename,
  SS_filename = SS_expected_filename,
  sigma_G_filename = sigma_G_filename,
  num_instruments = 77,
  num_observations = 108990,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 500)
)

# WARNING: The posterior computation takes a lot of time (days) and requires a lot of memory.
run_system_command(paste('./BayesMR', ini_filename))

Parkinson_BMI_posterior <- 
  read_PolyChord_samples(sprintf('chains/%s_temp_equal_weights.txt', fileroot))

# Save data files
save(Parkinson_BMI_posterior, BMI_genetic_associations,
     file = 'data/Parkinson_BMI_reproduced.RData')


# 6.4 Does coffee consumption influence smoking? - Figure 22 --------------

# Extract BayesMR input from summary statistics
n <- 2
J <- 8

EAF <- c(0.61, 0.89, 0.63, 0.29, 0.28, 0.81, 0.33, 0.45)
sigma_G <- sqrt(n * EAF * (1 - EAF))

# Ware et al. 2017 (https://onlinelibrary.wiley.com/doi/pdf/10.1111/add.13888), Table 3 (Page 1849) + Table S3
# Coffee consumption (coffee cups) and smoking heaviness (cigarettes per day) in UK Biobank

gamma_hat <- c(0.03, 0.03, 0.05, 0.06, 0.05, 0.03, 0.09, 0.03)
gamma_se <- c(0.01, 0.02, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01)
gamma_N <- 30062

Gamma_hat <- c(-0.214, 0.203, 0.016, -0.006, 0.015, 0.207, 0.076, -0.144)
Gamma_se <- c(0.132, 0.202, 0.128, 0.140, 0.137, 0.167, 0.135, 0.125)

Gamma_N <- obs_N <- 8072
beta_hat <- 0.45

N <- min(gamma_N, Gamma_N, obs_N)

# Writing the PolyChord configuration file

fileroot <- 'coffee_consumption_smoking'

ini_filename <- sprintf('ini/%s.ini', fileroot)
SS_expected_filename <- sprintf('data/%s_SS_expected.txt', fileroot)
SS_reverse_filename <- sprintf('data/%s_SS_reverse.txt', fileroot)
sigma_G_filename <- sprintf('data/%s_sigma_G.txt', fileroot)

SS_expected <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, obs_N, EAF)
SS_reverse <- SS_expected[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]

write.table(SS_expected, file = SS_expected_filename, row.names = FALSE, col.names = FALSE)
write.table(SS_reverse, file = SS_reverse_filename, row.names = FALSE, col.names = FALSE)
write.table(sqrt(n * EAF * (1 - EAF)), file = sigma_G_filename, row.names = FALSE, col.names = FALSE)

# Run in expected direction (coffee consumption -> smoking), arbitrarily chosen
write_BayesMR_configuration_file(
  config_filename = ini_filename,
  SS_filename = SS_expected_filename,
  sigma_G_filename = sigma_G_filename,
  num_instruments = J,
  num_observations = N,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

# WARNING: The posterior computation takes a lot of time (days) and requires a lot of memory.
run_system_command(paste('./BayesMR', ini_filename))
coffee_consumption_smoking_direction_log_evidence <- 
  read_PolyChord_log_evidence(sprintf('chains/%s.stats', fileroot))
coffee_consumption_smoking_direction_posterior <- 
  read_PolyChord_samples(sprintf('chains/%s_temp_equal_weights.txt', fileroot))

save(coffee_consumption_smoking_direction_log_evidence, coffee_consumption_smoking_direction_posterior,
     file = 'data/coffee_consumption_smoking_direction_reproduced.RData')

# Run in reverse direction (smoking -> coffee consumption)
write_BayesMR_configuration_file(
  config_filename = ini_filename,
  SS_filename = SS_reverse_filename,
  sigma_G_filename = sigma_G_filename,
  num_instruments = J,
  num_observations = N,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

# WARNING: The posterior computation takes a lot of time (days) and requires a lot of memory.
run_system_command(paste('./BayesMR', ini_filename))
smoking_coffee_consumption_direction_log_evidence <- 
  read_PolyChord_log_evidence(sprintf('chains/%s.stats', fileroot))
smoking_coffee_consumption_direction_posterior <- 
  read_PolyChord_samples(sprintf('chains/%s_temp_equal_weights.txt', fileroot))

save(smoking_coffee_consumption_direction_log_evidence, smoking_coffee_consumption_direction_posterior,
     file = 'data/smoking_coffee_consumption_direction_reproduced.RData')

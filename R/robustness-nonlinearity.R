# Linear case
set.seed(1918)
samples <- 1e4
MAF <- 0.3
n <- 2
G <- rbinom(samples, n, MAF)
U <- rnorm(samples)
X <- G + 0.1 * U + rnorm(samples)
Y <- 0.1 * G + X + 0.1 * U + rnorm(samples)

Z <- cbind(1, G, X, Y)
SSS <- t(Z) %*% Z / samples
RSS <- SSS[c(1, 2, 4, 3), c(1, 2, 4, 3)]
sigma_G <- sd(G) # sigma_G <- sqrt(MAF * (1 - MAF) * n)

file_wpwc_SSc <- 'data/weak_ply_weak_conf_A=Inf_SSc.txt'
file_wpwc_SSr <- 'data/weak_ply_weak_conf_A=Inf_SSr.txt'
file_wpwc_SiG <- 'data/weak_ply_weak_conf_A=Inf_SiG.txt'

write.table(SSS, file = file_wpwc_SSc, row.names = FALSE, col.names = FALSE)
write.table(RSS, file = file_wpwc_SSr, row.names = FALSE, col.names = FALSE)
write.table(sigma_G, file = file_wpwc_SiG, row.names = FALSE, col.names = FALSE)

# Introduce parametric nonlinearity

parametric_nonlinear <- function(x, A = 0.5) {
  A * tanh(x / A)
}


for (A in 2^seq(-1, 5)) {
  set.seed(1918)
  
  samples <- 1e4
  MAF <- 0.3
  n <- 2
  J <- 1 # number of instruments
  G <- rbinom(samples, n, MAF)
  U <- rnorm(samples)
  X <- G + 0.1 * U + rnorm(samples)
  Y <- 0.1 * G + parametric_nonlinear(X, A) + 0.1 * U + rnorm(samples)
  
  Z <- cbind(1, G, X, Y)
  SSS <- t(Z) %*% Z / samples
  RSS <- SSS[c(1, 2, 4, 3), c(1, 2, 4, 3)]
  sigma_G <- sd(G) # sigma_G <- sqrt(MAF * (1 - MAF) * n)

  # sigma_G stays the same
  filename_dir <- sprintf('data/weak_ply_weak_conf_A=%3.1f_SSc.txt', A)
  filename_rev <- sprintf('data/weak_ply_weak_conf_A=%3.1f_SSr.txt', A)
  
  write.table(SSS, file = filename_dir, row.names = FALSE, col.names = FALSE)
  write.table(RSS, file = filename_rev, row.names = FALSE, col.names = FALSE)
  
  
  write_model_configuration_file(sprintf('ini/nonlinear_A=%3.1f_model_dir.ini', A),
                                 SS_filename = filename_dir, 
                                 sigmaG_filename = file_wpwc_SiG,
                                 num_instruments = J,
                                 num_observations = samples)
  
  write_PolyChord_configuration_file(
    sprintf('ini/nonlinear_A=%3.1f_dir.ini', A),
    num_parameters = 2 * J + 7,
    num_live_points = 2000,
    num_repeats = 25,
    output_file_root = sprintf("nonlinear_A=%3.1f_dir", A)
  )
  
  system(paste('bin/polychord_MR', sprintf('ini/nonlinear_A=%3.1f_dir.ini', A), sprintf('ini/nonlinear_A=%3.1f_model_dir.ini', A)))
  
  write_model_configuration_file(sprintf('ini/nonlinear_A=%3.1f_model_rev.ini', A),
                                 SS_filename = filename_rev, 
                                 sigmaG_filename = file_wpwc_SiG,
                                 num_instruments = J,
                                 num_observations = samples)
  
  write_PolyChord_configuration_file(
    sprintf('ini/nonlinear_A=%3.1f_rev.ini', A),
    num_parameters = 2 * J + 7,
    num_live_points = 2000,
    num_repeats = 25,
    output_file_root = sprintf("nonlinear_A=%3.1f_rev", A)
  )    
  
  system(paste('bin/polychord_MR', sprintf('ini/nonlinear_A=%3.1f_rev.ini', A), sprintf('ini/nonlinear_A=%3.1f_model_rev.ini', A)))
}





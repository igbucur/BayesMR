source('R/prior_functions.R')

read_polychord_logevidence <- function(filename) {
  line <- readLines(filename, 9)[9]
  lev <-  as.numeric(substr(line, 17, 39))
  se <- as.numeric(substr(line, 47, 68))
  list(lev = lev, se = se,
       lower = lev - 2 * se, upper = lev + 2 * se)
}

derive_polychord_Bayes_factor <- function(filename_dir, filename_rev) {
  lev_dir <- read_polychord_logevidence(filename_dir)
  lev_rev <- read_polychord_logevidence(filename_rev)
  
  Bf <- exp(lev_dir$lev - lev_rev$lev)
  lower <- exp(lev_dir$lower - lev_rev$upper)
  upper <- exp(lev_dir$upper - lev_rev$lower)
  prob_dir <- Bf / (1 + Bf)
  prob_rev <- 1 - prob_dir
  
  list(
    Bf = Bf, lower = lower, upper = upper,
    prob_dir = prob_dir, prob_rev = prob_rev
  )
}

quick_derive_beta <- function(
  data, idx_sbeta = ncol(data) - 4, idx_sigma_X = ncol(data) - 1, idx_sigma_Y = ncol(data), 
  slab = 1, spike = 1e2
) {
  sigma_X <- qnorm(data[, idx_sigma_X] / 2 + 0.5, sd = 10)
  sigma_Y <- qnorm(data[, idx_sigma_Y] / 2 + 0.5, sd = 10)
  sbeta <- sapply(data[, idx_sbeta], quantile_spike_and_slab_2, slab = slab, spike = spike)
  
  sbeta * sigma_Y / sigma_X
}

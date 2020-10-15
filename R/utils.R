#' A more portable function for invoking a system command
#'
#' @param command String containing the system command to be invoked.
#' @param Windows_shell Full path to the Windows shell to be used for invoking 
#' the system command. MSYS2+MinGW has been tested and is recommended.
#' @param ... Further shell parameters.
#'
#' @return See \link[base]{system} for Linux and \link[base]{shell} for Windows.
#' @export
#'
#' @examples run_system_command("echo Hello!")
run_system_command <- function(command, Windows_shell = 'C:/msys64/msys2_shell.cmd -defterm -here -no-start -mingw64', ...) {
  
  if(.Platform$OS.type == "windows") {
    shell(shQuote(command), Windows_shell, flag = "-c", ...)
  } else {
    system(command)
  }
}


#' Function to generate data from a BayesMR generating model.
#'
#' @param N Integer number of observations.
#' @param theta Numeric expected allele frequencies for the genetic variables.
#' @param gamma Numeric direct causal effect from genetic variants to exposure.
#' @param alpha Numeric direct causal effect from genetic variants to outcome.
#' @param beta Numeric causal effect from exposure to outcome.
#' @param kappa_X Numeric confounding effect from confounder to exposure.
#' @param kappa_Y Numeric confounding effect from confounder to outcome.
#' @param sigma_X Numeric intrinsic standard deviation of exposure.
#' @param sigma_Y Numeric intrinsic standard deviation of outcome.
#' @param seed Integer random seed for reproducibility.
#' @param n Integer number of alleles (trials) for the binomial genetic variables.
#'
#' @return List containing generate data vector as well as the scatter matrix
#' of first-order and second-order statistics.
#' @export
#'
#' @examples
#' generate_data_BayesMR_model(N = 1000, theta = c(0.5, 0.3), gamma = c(1, 1), 
#' alpha = c(0, 0.05), beta = 1, kappa_X = 1, kappa_Y = 1, sigma_X = 1, sigma_Y = 1)
generate_data_BayesMR_model <- function(
  N, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, seed = NULL, n = 2
  ) {
  
  set.seed(seed) # set random seed
  
  G <- sapply(theta, function(t) stats::rbinom(N, n, t)) # genetic variant
  U <- stats::rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + stats::rnorm(N, sd = sigma_X) # exposure
  Y <- G %*% alpha + kappa_Y * U + beta * X + stats::rnorm(N, sd = sigma_Y) # outcome
  Z <- cbind(1, G, X, Y) # vector containing all variables
  
  list(data = Z, SS = t(Z) %*% Z / N)
}





#' Function for computing spike-and-slab quantiles.
#'
#' @param p Numeric probability value between 0 and 1.
#' @param w Numeric weight of spike-and-slab mixture between 0 (spike) and 1 (slab).
#' @param slab_precision Numeric precision of slab component.
#' @param spike_precision Numeric precision of spike component.
#'
#' @return Spike-and-slab quantile for given probability value.
#' @export
#'
#' @examples
#' quantile_spike_and_slab(0.9)
#' quantile_spike_and_slab(0.5, w = 0.3)
quantile_spike_and_slab <- function(p, w = 0.5, slab_precision = 1, spike_precision = 1e2) {
  
  qspike <- function(p) stats::qnorm(p, 0, sd = sqrt(1 / spike_precision))
  qslab <- function(p) stats::qnorm(p, 0, sd = sqrt(1 / slab_precision))
  pspike <- function(q) stats::pnorm(q, 0, sd = sqrt(1 / spike_precision))
  pslab <- function(q) stats::pnorm(q, 0, sd = sqrt(1 / slab_precision))
  
  
  left <- ifelse(p < 0.5, qslab(p), qspike(p))
  left_value <- w * pslab(left) + (1 - w) * pspike(left) - p
  
  right <- ifelse(p < 0.5, qspike(p), qslab(p))
  right_value <- w * pslab(right) + (1 - w) * pspike(right) - p
  
  middle <- (left + right) / 2
  mid_value <- w * pslab(middle) + (1 - w) * pspike(middle) - p
  
  
  while(abs(mid_value) > 1e-12) {
    if (left_value * mid_value < 0) {
      right <- middle
      right_value <- mid_value
    } else {
      left <- middle
      left_value <- mid_value
    }
    middle = (left + right) / 2
    mid_value = w * pslab(middle) + (1 - w) * pspike(middle) - p
  }
  middle
}

#' Vectorized version of function for computing spike-and-slab prior quantiles.
#'
#' @param p_vec Numeric vector of probability values between 0 and 1. 
#' @param x_vec Numeric vector of weights for spike-and-slab mixtures, all
#' between 0 (spike) and 1 (slab). Alternatively, a single value to be used for 
#' all p in p_vec.
#' @param slab_vec Numeric vector of slab component precisions. Alternatively,
#' a single value to be used for all p in p_vec.
#' @param spike_vec Numeric vector of spike component precisions. Alternatively,
#' a single value to be used for all p in p_vec.
#'
#' @return Vector of spike-and-slab quantile for given probability values.
#' @export
#'
#' @examples
#' quantile_spike_and_slab_vectorized(c(0.9, 0.5), c(0.5, 0.3))
quantile_spike_and_slab_vectorized <- function(p_vec, x_vec, slab_vec = 1, spike_vec = 1e2) {
  
  if (length(x_vec) == 1) x_vec <- rep(x_vec, length(p_vec))
  if (length(slab_vec) == 1) slab_vec <- rep(slab_vec, length(p_vec))
  if (length(spike_vec) == 1) spike_vec <- rep(spike_vec, length(p_vec))
  
  sapply(1:length(p_vec), function(i) {
    quantile_spike_and_slab(p_vec[i], x_vec[i], slab_vec[i], spike_vec[i])
  })
}

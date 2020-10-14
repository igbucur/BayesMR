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
#' @param N 
#' @param theta 
#' @param gamma 
#' @param alpha 
#' @param beta 
#' @param kappa_X 
#' @param kappa_Y 
#' @param sigma_X 
#' @param sigma_Y 
#' @param seed
#' @param n  
#'
#' @return
#' @export
#'
#' @examples
generate_data_BayesMR_model <- function(
  N, theta, gamma, alpha, beta, kappa_X, kappa_Y, sigma_X, sigma_Y, seed = NULL, n = 2
  ) {
  
  set.seed(seed)
  
  G <- sapply(theta, function(t) rbinom(N, n, t)) # genetic variant
  U <- rnorm(N) # confounder
  X <- G %*% gamma + kappa_X * U + rnorm(N, sd = sigma_X)
  Y <- G %*% alpha + kappa_Y * U + beta * X + rnorm(N, sd = sigma_Y)
  Z <- cbind(1, G, X, Y)
  
  list(data = Z, SS = t(Z) %*% Z / N)
}





# Get quantile for the spike and slab prior
quantile_spike_and_slab <- function(p, w = 0.5, slab = 1, spike = 1e2) {
  
  qspike <- function(p) stats::qnorm(p, 0, sd = sqrt(1 / spike))
  qslab <- function(p) stats::qnorm(p, 0, sd = sqrt(1 / slab))
  pspike <- function(q) stats::pnorm(q, 0, sd = sqrt(1 / spike))
  pslab <- function(q) stats::pnorm(q, 0, sd = sqrt(1 / slab))
  
  
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

quantile_spike_and_slab_vectorized <- function(p_vec, x_vec, slab_vec = 1, spike_vec = 1e2) {
  
  if (length(x_vec) == 1) x_vec <- rep(x_vec, length(p_vec))
  if (length(slab_vec) == 1) slab_vec <- rep(slab_vec, length(p_vec))
  if (length(spike_vec) == 1) spike_vec <- rep(spike_vec, length(p_vec))
  
  sapply(1:length(p_vec), function(i) {
    quantile_spike_and_slab(p_vec[i], x_vec[i], slab_vec[i], spike_vec[i])
  })
}

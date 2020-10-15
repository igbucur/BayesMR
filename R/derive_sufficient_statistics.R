

#' Function to extract first-order and second-order estimated moments from
#' summary data.
#'
#' @param J Integer number of instrumental variables.
#' @param beta_XG Numeric vector of estimated G -> X regression terms.
#' @param sigma_XG Numeric vector of estimated G -> X standard errors.
#' @param obs_XG Integer sample size in G -> X regression (GWA) study.
#' @param beta_YG Numeric vector of estimated G -> Y regression terms.
#' @param sigma_YG Numeric vector of estimated G -> Y standard errors.
#' @param obs_YG Integer sample size in G -> Y regression (GWA) study.
#' @param beta_YX Numeric vector of estimated X -> Y regression value.
#' @param obs_YX Integer sample size in X -> Y regression (observational) study.
#' @param EAF Numeric vector of effect allele frequencies
#' @param n Integer number of alleles (trials) for binomial genetic variant.
#'
#' @return Numeric matrix of size (J+3)x(J+3) containing the first-order and 
#' second-order moments of the (J+2) vector (G, X, Y).
#' @export
#'
#' @examples
#' MR_regression_coefficients_to_moments(1, 1, 1e-3, 1000, 1, 1e-3, 1000, 1, 1e-3)
MR_regression_coefficients_to_moments <- function(J, beta_XG, sigma_XG, obs_XG, beta_YG, sigma_YG, obs_YG, beta_YX, obs_YX, EAF = rep(0.5, J), n = 2) {

  N <- min(obs_YX, obs_XG, obs_YG)
  
  SSS <- matrix(0, J + 3, J + 3) # EV{[1 G X Y] [1 G X Y]^T}
  
  SSS[1, 1] <- 1
  SSS[1, 2:(J+1)] <- SSS[2:(J+1), 1] <- n * EAF
  SSS[2:(J+1), 2:(J+1)] <- n * n * EAF %*% t(EAF) + diag(n * EAF * (1 - EAF), J) # E{GG^T}
  
  SSS[1, J+2] <- SSS[J+2, 1] <- t(beta_XG) %*% SSS[1, 2:(J+1)]
  SSS[1, J+3] <- SSS[J+3, 1] <- t(beta_YG) %*% SSS[1, 2:(J+1)]
  
  SSS[J+2, 2:(J+1)] <- SSS[2:(J+1), J+2] <- SSS[2:(J+1), 2:(J+1)] %*% beta_XG
  SSS[J+3, 2:(J+1)] <- SSS[2:(J+1), J+3] <- SSS[2:(J+1), 2:(J+1)] %*% beta_YG
  
  SSS[J+2, J+2] <- mean((beta_XG^2 + sigma_XG^2 * obs_XG) * n * EAF * (1 - EAF)) + (n * t(beta_XG) %*% EAF)^2
  SSS[J+2, J+3] <- SSS[J+3, J+2] <- mean((beta_XG^2 + sigma_XG^2 * obs_XG) * n * EAF * (1 - EAF)) * beta_YX + SSS[1, J+2] * SSS[J+3, 1]
  SSS[J+3, J+3] <- mean((beta_YG^2 + sigma_YG^2 * obs_YG) * n * EAF * (1 - EAF)) + (n * t(beta_YG) %*% EAF)^2
  
  SSS
}

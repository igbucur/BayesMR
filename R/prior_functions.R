
#' @title horseshoe_prior
#' @description Implementation of Horseshoe Prior for a single causal effect
#' @param b causal effect variable
#' @return prior pdf value at input
horseshoe_prior <- function(b) {
  # use approximation
  log(log(1 + 2 / (b * b))) # - 0.5 * log(2 * pi * pi * pi)
}

log_spike_and_slab <- function(b, w, k_slab, k_spike) {
  
  sd_slab <- sqrt(1 / k_slab)
  sd_spike <- sqrt(1 / k_spike)

  log(w * dnorm(b, 0, sd_slab) + (1 - w) * dnorm(b, 0, sd_spike))
}


#' Robust version of log spike-and-slab prior
#'
#' @param b parameter for which we compute the prior
#' @param w weight of slab vs spike
#' @param k_slab precision of spike 
#' @param k_spike precision of slab
#'
#' @return
#' @export
#'
#' @examples
log_spike_and_slab_robust <- function(b, w, k_slab, k_spike) {

  sd_slab <- sqrt(1 / k_slab)
  sd_spike <- sqrt(1 / k_spike)
  
  -0.5 * (k_slab * b * b + log(2 * pi)) + 
    log(w * sqrt(k_slab) + (1 - w) * sqrt(k_spike) * exp(- 0.5 * (k_spike - k_slab) * b * b))
}


#' Robust version of log spike-and-slab prior for matrices
#'
#' @param B lower triangular matrix of parameters for which to compute log spike-and-slab priors
#' @param w weight of slab vs spike
#' @param k_slab precision of spike 
#' @param k_spike precision of slab
#'
#' @return sum of log spike-and-slab priors for every parameter in matrix B
#' @export
#'
#' @examples
log_spike_and_slab_matrix_robust <- function(B, w, k_slab, k_spike) {
  
  b <- B[lower.tri(B)]
  
  sd_slab <- sqrt(1 / k_slab)
  sd_spike <- sqrt(1 / k_spike)
  
  sum(log_spike_and_slab_robust(b, w, k_slab, k_spike))
}


log_prior_variance <- function(V) {
  v <- diag(V)
  
  sum(-log(1 + v * v))
}

log_prior_confounding <- function(C, k) {
  p <- nrow(C)
  q <- ncol(C)
  
  stopifnot(q == p * (p - 1) / 2)

  
  # 2 x q matrix
  edges <- combn(p, 2)
  
  sum(dnorm(C[cbind(edges[1,], 1:q)], sd = sqrt(1 / k), log = TRUE)) +
    sum(dnorm(C[cbind(edges[2,], 1:q)], sd = sqrt(1 / k), log = TRUE))
}

logPrior <- function(B, V, beta = 0.5, k1 = 1, k2 = 1e3 * k1) {
  b <- B[lower.tri(B)]
  v <- diag(V)
  
  sd1 <- sqrt(1 / k1)
  sd2 <- sqrt(1 / k2)
  
  sum(c(log(beta * dnorm(b, 0, sd1) + (1 - beta) * dnorm(b, 0, sd2)), -log(v)))
}

logPriorHalfCauchy <- function(B, V, beta = 0.5, k1 = 1, k2 = 1e3 * k1) {
  b <- B[lower.tri(B)]
  v <- diag(V)
  
  sd1 <- sqrt(1 / k1)
  sd2 <- sqrt(1 / k2)
  
  sum(c(log(beta * dnorm(b, 0, sd1) + (1 - beta) * dnorm(b, 0, sd2)), -log(1 + v * v)))
}

# \alpha * \mathcal{N}(0, \sigma^2) + (1 - \alpha) * \mathcal{N}(0, c\sigma^2)
logS3Prior <- function(alpha, vtail, varrow, beta = 0.5, k1 = 1, k2 = 1e2) {

  # r <- sqrt(vtail / varrow)
  sd1 <- sqrt(varrow / vtail / k1)
  sd2 <- sqrt(varrow / vtail / k2)
  log(beta * dnorm(alpha, 0, sd1) + (1 - beta) * dnorm(alpha, 0, sd2))
  # dnorm(alpha, 0, sdspike * r, log = TRUE)
}


#' Get quantile for the spike and slab prior
quantile_spike_and_slab <- function(p, beta = 0.5, slab = 1, spike = 1e4) {
  
  if (slab == spike) return (qnorm(p, 0, sd = sqrt(1 / slab)))
  if (abs(beta - 1) < 1e-12) return (qnorm(p, 0, sd = sqrt(1 / slab)))
  if (abs(beta) < 1e-12) return (qnorm(p, 0, sd = sqrt(1 / spike)))
  
  bound1 <- qnorm(p, 0, sd = sqrt(1 / spike))
  bound2 <- qnorm(p, 0, sd = sqrt(1 / slab))
  
  if (bound1 == bound2) return (bound1)
  
  f <- function(q) {
    beta * pnorm(q, 0, sd = sqrt(1 / slab)) + (1 - beta) * pnorm(q, 0, sd = sqrt(1 / spike)) - p
  }
  
  uniroot(f, c(min(bound1, bound2), max(bound1, bound2)), tol = 1e-12)$root
}

quantile_spike_and_slab_2 <- function(p, beta = 0.5, slab = 1, spike = 1e4) {
  
  qspike <- function(p) qnorm(p, 0, sd = sqrt(1 / spike))
  qslab <- function(p) qnorm(p, 0, sd = sqrt(1 / slab))
  pspike <- function(q) pnorm(q, 0, sd = sqrt(1 / spike))
  pslab <- function(q) pnorm(q, 0, sd = sqrt(1 / slab))
  
  
  left <- ifelse(p < 0.5, qslab(p), qspike(p))
  left_value <- beta * pslab(left) + (1 - beta) * pspike(left) - p
  
  right <- ifelse(p < 0.5, qspike(p), qslab(p))
  right_value <- beta * pslab(right) + (1 - beta) * pspike(right) - p
  
  middle <- (left + right) / 2
  mid_value <- beta * pslab(middle) + (1 - beta) * pspike(middle) - p
  
  
  while(abs(mid_value) > 1e-12) {
    if (left_value * mid_value < 0) {
      right <- middle
      right_value <- mid_value
    } else {
      left <- middle
      left_value <- mid_value
    }
    middle = (left + right) / 2
    mid_value = beta * pslab(middle) + (1 - beta) * pspike(middle) - p
  }
  middle
}

#' Get quantile for b43 given b42: (b42, b43) are distributed according to multivariate spike-and-slab prior
quantile_spike_and_slab_b43 <- function(p, b42, beta = 0.5, slab = 1, spike = 1e4) {
  
  beta <- beta * dnorm(b42, 0, sd = sqrt(1 / slab)) / (beta * dnorm(b42, 0, sd = sqrt(1 / slab)) + (1 - beta) * dnorm(b42, 0, sd = sqrt(1 / spike)))
  
  quantile_spike_and_slab(p, beta, slab, spike)
}


log_prior_original_scale <- function(b, v_tail, v_arrow, w_slab = 0.5, v_slab = 1, v_spike = 1e-2) {
  
  logS3Prior(b, v_tail, v_arrow, w_slab, 1 / v_slab, 1 / v_spike)
}

log_prior_after_rescaling <- function(b, w_slab = 0.5, v_slab = 1, v_spike = 1e-2) {
  sd_slab <- sqrt(v_slab)
  sd_spike <- sqrt(v_spike)
  
  log(w_slab * dnorm(b, sd = sd_slab) + (1 - w_slab) * dnorm(b, sd = sd_spike))
}

# 
# logPriorGradient <- function(alpha, beta, sdspike, sdslab) {
#   - alpha * (beta * dnorm(alpha, 0, sdspike) / (sdspike * sdspike) + (1 - beta) * dnorm(alpha, 0, sdslab) / (sdslab * sdslab)) / 
#     priorDensityMixtureGaussian(alpha, beta, sdspike, sdslab)
# }
# 
# logPriorHessian <- function(alpha, beta, sdspike, sdslab) {
#   (-(beta * dnorm(alpha, 0, sdspike) / (sdspike * sdspike) + (1 - beta) * dnorm(alpha, 0, sdslab) / (sdslab * sdslab)) * priorDensityMixtureGaussian(alpha, beta, sdspike, sdslab) +
#      beta * (1 - beta) * alpha * alpha * dnorm(alpha, 0, sdspike) * dnorm(alpha, 0, sdslab) * (1 / (sdspike * sdspike) -  1 / (sdslab * sdslab))^2) / 
#     (priorDensityMixtureGaussian(alpha, beta, sdspike, sdslab) * priorDensityMixtureGaussian(alpha, beta, sdspike, sdslab))
# }

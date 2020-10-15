source('R/utils.R')

#' Read data.frame from a PolyChord output file and name the columns appropriately.
#'
#' @param filename Should be of the form {root}.txt or {root}_equal_weights.txt (see PolyChord documentation)
#' @param slab_precision Numeric precision of slab component in the BayesMR priors.
#' @param spike_precision Numeric precision of spike component in the BayesMR priors.
#' @param quantile_transform Logical flag (true by default) indicating whether
#' the samples read from the PolyChord output should be transformed back from the
#' uniform hypercube to the parameters' distributions via the quantile function.
#'
#' @return A data frame containing the PolyChord samples for the BayesMR model variables.
#' @export
read_PolyChord_samples <- function(filename, slab_precision = 1, spike_precision = 100,
                                   quantile_transform = TRUE) {
  
  data_frame <- utils::read.table(filename)
  J <- (ncol(data_frame) - 9) / 2 # number of instruments
  names(data_frame) <- c(
    'w', 'lp__', 'wgamma', 'walpha',
    paste0(c(outer(c('sgamma[', 'salpha['), 1:J, paste0)), ']'),
    'sbeta', 'skappa_X', 'skappa_Y', 'sigma_X', 'sigma_Y'
  )
  
  if (quantile_transform) {
    for (j in 1:J) {
      sgamma <- paste0('sgamma[', j, ']')
      salpha <- paste0('salpha[', j, ']')
      data_frame[, sgamma] <- quantile_spike_and_slab_vectorized(
        data_frame[, sgamma], data_frame[, 'wgamma'], slab_precision, spike_precision
      )
      data_frame[, salpha] <- quantile_spike_and_slab_vectorized(
        data_frame[, salpha], data_frame[, 'walpha'], slab_precision, spike_precision
      )
    }
    
    data_frame[, 'sbeta'] <- quantile_spike_and_slab_vectorized(
      data_frame[, 'sbeta'], 0.5, slab_precision, spike_precision
    )
    
    data_frame[, 'skappa_X'] <- quantile_spike_and_slab_vectorized(
      data_frame[, 'skappa_X'] / 2 + 0.5, 0.5, slab_precision, spike_precision
    )
    
    data_frame[, 'skappa_Y'] <- quantile_spike_and_slab_vectorized(
      data_frame[, 'skappa_Y'], 0.5, slab_precision, spike_precision
    )
    
    data_frame[, 'sigma_X'] <- stats::qnorm(data_frame[, 'sigma_X'] / 2 + 0.5, sd = 10)
    data_frame[, 'sigma_Y'] <- stats::qnorm(data_frame[, 'sigma_Y'] / 2 + 0.5, sd = 10)
  }
  
  data_frame
}

#' Function for reading the log-evidence from PolyChord output.
#'
#' @param filename Name of the PolyChord file
#'
#' @return List containing estimated log-evidence, standard error, as well
#' as derived lower and upper 95% CI error bars.
#' @export
read_PolyChord_log_evidence <- function(filename) {
  line <- readLines(filename, 9)[9]
  lev <-  as.numeric(substr(line, 17, 39))
  se <- as.numeric(substr(line, 47, 68))
  list(lev = lev, se = se,
       lower = lev - 2 * se, upper = lev + 2 * se)
}

#' Function for comparing the log-evidence from two BayesMR runs corresponding
#' to both directions of the causal link and deriving the Bayes factor.
#'
#' @param filename_expected Name of PolyChord output file containing log-evidence
#' for the BayesMR run in the expected causal link direction.
#' @param filename_reverse Name of PolyChord output file containing log-evidence
#' for the BayesMR run in the reverse causal link direction.
#'
#' @return List containing Bayes factor together with its lower and upper bar,
#' as well as the probability of the expected and reverse direction.
#' @export
derive_polychord_Bayes_factor <- function(filename_expected, filename_reverse) {
  
  log_evidence_expected <- read_PolyChord_log_evidence(filename_expected)
  log_evidence_reverse <- read_PolyChord_log_evidence(filename_reverse)
  
  Bf <- exp(log_evidence_expected$lev - log_evidence_reverse$lev)
  lower <- exp(log_evidence_expected$lower - log_evidence_reverse$upper)
  upper <- exp(log_evidence_expected$upper - log_evidence_reverse$lower)
  prob_dir <- Bf / (1 + Bf)
  prob_rev <- 1 - prob_dir
  
  list(
    Bf = Bf, lower = lower, upper = upper,
    prob_dir = prob_dir, prob_rev = prob_rev
  )
}

#' Function for extracting the causal effect from exposure to outcome from a
#' data frame of PolyChord samples.
#'
#' @param data Data frame containing PolyChord samples, typically read using read.table.
#' @param idx_sbeta Integer index of sbeta variable.
#' @param idx_sigma_X Integer index of sigma_X variable.
#' @param idx_sigma_Y Integer index of sigma_Y variable.
#' @param slab_precision Numeric precision of slab component.
#' @param spike_precision Numeric precision of spike component.
#'
#' @return Numeric vector containing samples of beta (the causal effect from
#' exposure to outcome) estimates in the original scale.
#' @export
extract_beta_PolyChord_samples <- function(
  data, idx_sbeta = ncol(data) - 4, idx_sigma_X = ncol(data) - 1, idx_sigma_Y = ncol(data), 
  slab_precision = 1, spike_precision = 1e2
) {
  sigma_X <- stats::qnorm(data[, idx_sigma_X] / 2 + 0.5, sd = 10)
  sigma_Y <- stats::qnorm(data[, idx_sigma_Y] / 2 + 0.5, sd = 10)
  sbeta <- sapply(data[, idx_sbeta], quantile_spike_and_slab, 
                  slab_precision = slab_precision, spike_precision = spike_precision)
  
  sbeta * sigma_Y / sigma_X
}

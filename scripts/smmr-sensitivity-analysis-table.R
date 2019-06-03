


nonlinear_evidence <- list()

for (A in 2^seq(-1, 5)) {
  nonlinear_evidence[[as.character(A)]] <- derive_polychord_Bayes_factor(
    sprintf('chains/nonlinear_A=%3.1f_dir.stats', A),
    sprintf('chains/nonlinear_A=%3.1f_rev.stats', A)
  )
}

samples_nonlinear <- list()
for (A in 2^seq(-1, 5)) {
  samples_nonlinear[['dir']][[as.character(A)]] <- read.table(sprintf('chains/nonlinear_A=%3.1f_dir_equal_weights.txt', A))
  samples_nonlinear[['rev']][[as.character(A)]] <- read.table(sprintf('chains/nonlinear_A=%3.1f_rev_equal_weights.txt', A))
}

beta_samples_nonlinear_dir <- lapply(samples_nonlinear$dir, quick_derive_beta)
#beta_samples_nonlinear_rev <- lapply(samples_nonlinear$rev, quick_derive_beta)

plot(sapply(nonlinear_evidence, '[[', 'prob_dir') ~ seq(-1, 5), type = 'l', ylim = c(0, 1),
     main = "Probability of expected direction for non-linear exposure-outcome relation",
     ylab = "Posterior probability estimate", xlab = expression(log[2](A)))

apply(prob_ci_nonlinear <- sapply(nonlinear_evidence, function(ev) {
  prob_lwr <- ev$lower / (1 + ev$lower)
  prob_upr <- ev$upper / (1 + ev$upper)
  
  c(prob_lwr, prob_upr)
}), 1, lines, x = seq(-1, 5), lty = 2)


nongaussian_evidence <- list()

for (A in 2^seq(0, 6)) {
  nongaussian_evidence[[as.character(A)]] <- derive_polychord_Bayes_factor(
    sprintf('chains/nongaussian_A=%3.1f_dir.stats', A),
    sprintf('chains/nongaussian_A=%3.1f_rev.stats', A)
  )
}

plot(sapply(nongaussian_evidence, '[[', 'prob_dir') ~ seq(0, 6), type = 'l', ylim = c(0, 1),
     main = "Probability of expected direction in the presence of non-Gaussian noise",
     ylab = "Posterior probability estimate", xlab = expression(log[2](nu)))

apply(prob_ci_nongaussian <- sapply(nongaussian_evidence, function(ev) {
  prob_lwr <- ev$lower / (1 + ev$lower)
  prob_upr <- ev$upper / (1 + ev$upper)
  
  c(prob_lwr, prob_upr)
}), 1, lines, x = seq(0, 6), lty = 2)

samples_nongaussian <- list()
for (A in 2^seq(0, 6)) {
  samples_nongaussian[['dir']][[as.character(A)]] <- read.table(sprintf('chains/nongaussian_A=%3.1f_dir_equal_weights.txt', A))
  samples_nongaussian[['rev']][[as.character(A)]] <- read.table(sprintf('chains/nongaussian_A=%3.1f_rev_equal_weights.txt', A))
}

beta_samples_nongaussian_dir <- lapply(samples_nongaussian$dir, quick_derive_beta)
#beta_samples_nongaussian_rev <- lapply(samples_nongaussian$rev, quick_derive_beta)

library(tidyverse)

beta_densities_non_gaussian_dir <- lapply(beta_samples_nongaussian_dir, density)

plot_limits <- sapply(beta_densities_non_gaussian_dir, function(d) c(min_x = min(d$x), max_x = max(d$x), max_y = max(d$y)))

plot(beta_densities_non_gaussian_dir$`1`, ylim = c(0, max(plot_limits['max_y', ])),
     xlim = c(min(plot_limits['min_x', ]), max(plot_limits['max_x', ])))
lines(beta_densities_non_gaussian_dir$`2`)
lines(beta_densities_non_gaussian_dir$`4`)
lines(beta_densities_non_gaussian_dir$`8`)
lines(beta_densities_non_gaussian_dir$`16`)
lines(beta_densities_non_gaussian_dir$`32`)

library(Hmisc)

#model_prob_nonlinear <- sapply(nonlinear_evidence, '[[', 'prob_dir')
model_prob_nonlinear <- apply(prob_ci_nonlinear, 2, function(col) paste0('[', round(col[1], 3), ', ', round(col[2], 3), ']'))
beta_summary_nonlinear <- round(t(sapply(beta_samples_nonlinear_dir, summary)[2:5, ]), 3)
table_data_nonlinear <- cbind(as.numeric(rownames(beta_summary_nonlinear)), beta_summary_nonlinear, model_prob_nonlinear)
latex(table_data_nonlinear, rowname = NULL, file = "", col.just = rep("c", 6), 
      # cgroup = c("a", "\\hat{\\ce}}", "b"), n.cgroup = c(1, 4, 1),
      colheads = c("A", "Q1", "Median", "Mean", "Q3", "$\\hat{p}(\\model_{X \\to Y} \\given \\data)$"), digits = 3)


#model_prob_nongaussian <- sapply(nongaussian_evidence, '[[', 'prob_dir')
model_prob_nongaussian <- apply(prob_ci_nongaussian, 2, function(col) paste0('[', round(col[1], 3), ', ', round(col[2], 3), ']'))
beta_summary_nonlinear <- round(t(sapply(beta_samples_nonlinear_dir, summary)[2:5, ]), 3)
beta_summary_nongaussian <- round(t(sapply(beta_samples_nongaussian_dir, summary)[2:5, ]), 3)
table_data_nongaussian <- cbind(as.numeric(rownames(beta_summary_nongaussian)), beta_summary_nongaussian, model_prob_nongaussian)
latex(table_data_nongaussian, rowname = NULL, file = "", col.just = rep("c", 6),
      colheads = c("\\nu", "Q1", "Median", "Mean", "Q3", "$\\hat{p}(\\model_{X \\to Y} \\given \\data)$"), digits = 3)

# 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(bayesplot)
library(knitr)
library(kableExtra)

source('R/PolyChord_interface.R') # for quick_derive_beta

figures_dir <- "figures/"
if (!dir.exists(figures_dir)) dir.create(figures_dir)


# 4.4 - Example: Near-conditional independence (near-LCD) - Figure 14 -----

load('data/near_LCD_example.RData')

model_averaging_plot <- function(example_data, expected_direction = TRUE, xlim = c(-0.25, 2)) {
  
  if (expected_direction) {
    beta_expected <- extract_beta_PolyChord_samples(example_data$samples_expected, spike = 100)
    fig <- ggplot(as.data.frame(beta_expected), aes(x = beta_expected)) +
      geom_density(aes_q(y = bquote(..density.. * .(example_data$evidence_comparison$probability_expected))), fill = 'gray') +
      theme_tufte(base_size = 30) + xlim(xlim) +
      geom_vline(xintercept = 0, lty = 'dashed') +
      xlab(expression(beta[X %->% Y])) +
      ylab("Posterior density")
    fig <- fig + annotate("text", x = xlim[2] / 12, y = 0.9 * max(ggplot_build(fig)$data[[1]]$y), 
               label = sprintf("%.2f%%", round(example_data$evidence_comparison$probability_reverse * 100, 2)), size = 8)
    
  } else {
    beta_reverse <- extract_beta_PolyChord_samples(example_data$samples_reverse, spike = 100)
    fig <- ggplot(as.data.frame(beta_reverse), aes(x = beta_reverse)) +
      geom_density(aes_q(y = bquote(..density.. * .(example_data$evidence_comparison$probability_reverse))), fill = 'gray') +
      theme_tufte(base_size = 30) + xlim(xlim) +
      geom_vline(xintercept = 0, lty = 'dashed') +
      annotate("text", x = annot_x, y = annot_y, label = sprintf("%.2f%%", round(example_data$evidence_comparison$probability_expected * 100, 2))) +
      xlab(expression(beta[Y %->% X])) +
      ylab("Posterior density")
    fig <- fig + annotate("text", x = xlim[2] / 12, y = 0.9 * max(ggplot_build(fig)$data[[1]]$y), 
               label = sprintf("%.2f%%", round(evidence_comparison$probability_expected * 100, 2)), size = 8)
  }
  
  fig
}

ggsave(paste0(figures_dir, "Figure_14.pdf"), 
       model_averaging_plot(near_LCD_example), 
       width = 10, height = 5)


# 5.1 - Estimation Robustness - Figure 15 ---------------------------------

load('data/estimation_robustness.RData')

estimation_robustness_plot <- mcmc_areas(estimation_robustness) +
  theme_tufte() + coord_cartesian(xlim = c(-0.5, 1.5)) +
  theme(text = element_text(size = 30), axis.title.x = element_text(margin = margin(t = 10))) +
  geom_vline(xintercept = 1) +
  xlab(expression(beta)) +
  ylab("Percentage of valid instruments")

ggsave(paste0(figures_dir, "Figure_15.pdf"), estimation_robustness_plot, 
       width = 10, height = 7)


# 5.2 - Sensitivity to Prior Hyperparameters - Figure 16 ------------------

load('data/sensitivity_to_hyperparameters.RData')

sensitivity_to_hyperparameters_plot <- 
  ggplot(sensitivity_to_hyperparameters, aes(x = spike, y = beta)) +
  geom_violin(fill = 'gray', colour = 'darkblue') +
  ylab(expression(beta))+
  xlab(expression(lambda)) +
  coord_flip(ylim = c(-0.5, 1.5)) +
  theme_tufte(base_size = 30) +
  theme(aspect.ratio = 2) +
  theme(axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = sapply(
    paste0("10^", -as.integer(levels(sensitivity_to_hyperparameters$spike))), 
    function(t) parse(text = t), USE.NAMES = FALSE)) +
  scale_linetype_manual(name = "", values = c("WPP bounds" = "dashed")) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  geom_hline(yintercept = 1, lty = 2)

ggsave(paste0(figures_dir, "Figure_16.pdf"), sensitivity_to_hyperparameters_plot, 
       height = 10, width = 6)


# 5.3 - Model Averaging Robustness - Table 1 ------------------------------

load('data/model_averaging_robustness.RData')

# Read data, combute beta estimates and their standard errors
tibble_bidirectional_MR <- model_averaging_robustness %>%
  mutate(expected_beta_hat = expected_Gamma_hat / expected_gamma_hat,
         expected_beta_std_err = expected_Gamma_std_err / expected_gamma_hat,
         reverse_beta_hat = reverse_Gamma_hat / reverse_gamma_hat,
         reverse_beta_std_err = reverse_Gamma_std_err / reverse_gamma_hat) %>%
  dplyr::select(delta, expected_beta_hat, expected_beta_std_err,
                reverse_beta_hat, reverse_beta_std_err, probability_expected)
  
knitr::kable(
  tibble_bidirectional_MR, digits = 3, format = "latex", booktabs = TRUE, escape = FALSE, 
  col.names = c("$\\delta$", "$\\hat{\\beta}^{IV}_{X \\to Y}$", "$\\hat{\\sigma}^{IV}_{X \\to Y}$",
    "$\\hat{\\beta}^{IV}_{Y \\to X}$", "$\\hat{\\sigma}^{IV}_{Y \\to X}$", 
    "$\\hat{p}(\\mathcal{M}_{X \\to Y} | \\textbf{D})$")) %>%
  kableExtra::kable_styling() %>%
  kableExtra::save_kable(paste0(figures_dir, "Table_1.pdf"))


# 5.4 - Sensitivity to Nonlinearity and Nonnormality - Tables 2 and 3 -----

load('data/sensitivity_to_parametric_assumptions.RData')

# Collect data in a tibble data frame
tibble_nonlinear <- tibble::as_tibble(t(sapply(sensitivity_to_nonlinearity, function(results) {
  c(round(summary(results$expected$beta_posterior)[2:5], 3),
    prob =paste0('[', round(results$direction_evidence$probability_expected_lower_error_bar, 3),
           ', ', round(results$direction_evidence$probability_expected_upper_error_bar, 3), ']'))
})), rownames = NA) %>%
  rownames_to_column() %>%
  mutate(rowname = gsub("A=", "", rowname))

# Create LaTex table and save it in PDF format
knitr::kable(
  tibble_nonlinear, format = "latex", row.names = FALSE, align = "c", booktabs = TRUE, digits = 3,
  col.names = c("A", "Q1", "Median", "Mean", "Q3", "$\\hat{p}(\\mathcal{M}_{X \\to Y} | \\textbf{D})$"), escape = FALSE
) %>% kableExtra::kable_styling() %>%
  kableExtra::add_header_above(header = c(" ", "$\\\\hat{\\\\beta}_{X \\\\to Y}$" = 2), 
                               escape = FALSE, align = "l") %>%
  kableExtra::save_kable(paste0(figures_dir, "Table_2.pdf"))

# Collect data in a tibble data frame
tibble_nonnormal <- tibble::as_tibble(t(sapply(sensitivity_to_nonnormality, function(results) {
  c(round(summary(results$expected$beta_posterior)[2:5], 3),
    prob =paste0('[', round(results$direction_evidence$probability_expected_lower_error_bar, 3),
                 ', ', round(results$direction_evidence$probability_expected_upper_error_bar, 3), ']'))
})), rownames = NA) %>%
  rownames_to_column() %>%
  mutate(rowname = gsub("nu=", "", rowname))

# Create LaTex table and save it in PDF format
knitr::kable(
  tibble_nonnormal, format = "latex", row.names = FALSE, align = "c", booktabs = TRUE, digits = 3,
  col.names = c("\\nu", "Q1", "Median", "Mean", "Q3", "$\\hat{p}(\\mathcal{M}_{X \\to Y} | \\textbf{D})$"), escape = FALSE
) %>% kableExtra::kable_styling() %>%
  kableExtra::add_header_above(header = c(" ", "$\\\\hat{\\\\beta}_{X \\\\to Y}$" = 2), 
                               escape = FALSE, align = "l") %>%
  kableExtra::save_kable(paste0(figures_dir, "Table_3.pdf"))

# 6.2 - Effect of birth weight on fasting glucose - Figures 17 and 18 -----

load('data/birth_weight_fasting_glucose.RData')

plot_birth_weight_fasting_glucose_beta <- function(samples) {
  mcmc_areas(samples, pars = "beta") +
  theme_tufte() + coord_cartesian(xlim = c(-0.4, 0.2)) + 
  theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.155) + geom_vline(xintercept = c(-0.223, -0.088), linetype = "dashed") + 
  xlab(expression(beta)) + ylab("Posterior density") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -4.3, vjust = -3.1, parse = T, size = 10)
}

# Figure 17 - BayesMR with default settings
ggsave(paste0(figures_dir, "Figure_17.pdf"),
       plot_birth_weight_fasting_glucose_beta(birth_weight_fasting_glucose_default),
       width = 10, height = 5)

# Figure 18 - BayesMR incorporating IV assumptions
ggsave(paste0(figures_dir, "Figure_18.pdf"), 
       plot_birth_weight_fasting_glucose_beta(birth_weight_fasting_glucose_IV), 
       width = 10, height = 5)



# 6.3 Effect of BMI on the risk of PD - Figures 19-21 ---------------------

load('data/Parkinson_BMI.RData')
genetic_associations <- read.csv('inst/extdata/Parkinson_BMI_genetic_associations.csv')

Parkinson_BMI_beta_plot <- mcmc_areas(Parkinson_BMI_posterior, pars = "beta") +
  ggthemes::theme_tufte() + 
  theme(text = element_text(size = 30), 
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.1946717) + 
  geom_vline(xintercept = c(-0.3682688, -0.02107458), linetype = "dashed") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", 
           hjust = -5.2, vjust = -10, parse = T, size = 8) +
  scale_y_discrete(labels = "") + ylab("Posterior density") +
  xlab(expression(paste("Effect of BMI on the risk of PD (", beta, ")")))

ggsave(paste0(figures_dir, "Figure_19.pdf"), Parkinson_BMI_beta_plot, width = 8, height = 8)

Parkinson_BMI_outliers_plot <- 
  ggplot(data.frame(genetic_associations), aes(x = gamma_hat, y = Gamma_hat)) + 
  geom_point(color = c('red', 'red', rep('black', nrow(genetic_associations) - 2)), 
             shape = c(2, 2, rep(1, nrow(genetic_associations) - 2)), size = 3) + 
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col = 'red', linetype = 2) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col = 'black',
              data = data.frame(genetic_associations)[-c(1, 2),]) +
  ggthemes::theme_tufte() + 
  theme(text = element_text(size = 25), 
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.title.y = element_text(margin = margin(r = 10))) +
  xlab("Genetic associations with the exposure") +
  ylab("Genetic associations with the outcome")

ggsave(paste0(figures_dir, "Figure_20.pdf"), Parkinson_BMI_outliers_plot, fonts = "serif", width = 10, height = 7)

Parkinson_BMI_alpha_plot <- mcmc_areas(Parkinson_BMI_posterior, pars = c("alpha[1]", "alpha[2]")) +
  theme(text = element_text(size = 30), axis.title.x = element_text(margin = margin(t = 30))) +
  ggthemes::theme_tufte() + 
  scale_y_discrete(labels = c("rs17001654", "rs13107325")) +
  theme(text = element_text(size = 25), 
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.title.y = element_text(margin = margin(r = 10))) +
  xlab(expression(paste("Pleiotropic effect on risk of PD (", alpha, ")")))

ggsave(paste0(figures_dir, "Figure_21.pdf"), Parkinson_BMI_alpha_plot, width = 10, height = 8)


# 6.4 Does coffee consumption influence smoking? - Figure 22 --------------

load('data/coffee_consumption_smoking_direction.RData')
load('data/smoking_coffee_consumption_direction.RData')

# evidence when spike = 1e-4, obtained with PolyChord 2000 live points
log_evidence_coffee_consumption_smoking_direction <- -0.422878161448718E+005
log_evidence_coffee_consumption_smoking_direction_standard_error <- 0.140144123185408E+000
log_evidence_smoking_coffee_consumption_direction <- -0.422876287259804E+005
log_evidence_smoking_coffee_consumption_direction_standard_error <- 0.133930400352899E+000

# Estimated probability of expected and reverse direction (in that order)
probability_direction <- 
  c(1, exp(log_evidence_smoking_coffee_consumption_direction - log_evidence_coffee_consumption_smoking_direction)) / 
  (1 + exp(log_evidence_smoking_coffee_consumption_direction - log_evidence_coffee_consumption_smoking_direction))

coffee_consumption_smoking_direction_plot <- 
  ggplot(coffee_consumption_smoking_direction_posterior, aes(x = beta)) +
  geom_density(aes_q(y = bquote(..density.. * .(probability_direction[1]))), fill = 'gray') +
  theme_tufte(base_size = 30) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  annotate("text", x = -0.25, y = 1.5, label = sprintf("%.2f%%", round(probability_direction[2] * 100, 2)), size = 8) +
  xlim(c(-1, 2)) +
  xlab("Causal effect of coffee consumption on smoking") +
  ylab("Posterior density")

ggsave(paste0(figures_dir, "Figure_22a.pdf"), 
       coffee_consumption_smoking_direction_plot, 
       width = 10, height = 7)

smoking_coffee_consumption_direction_plot <- 
  ggplot(smoking_coffee_consumption_direction_posterior, aes(x = beta)) +
  geom_density(aes_q(y = bquote(..density.. * .(probability_direction[2]))), fill = 'gray') +
  theme_tufte(base_size = 30) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  annotate("text", x = -0.03, y = 60, label = sprintf("%.2f%%", round(probability_direction[1] * 100, 2)), size = 8) +
  xlim(c(-0.2, 0.2)) +
  xlab("Causal effect of smoking on coffee consumption") +
  ylab("Posterior density")

ggsave(paste0(figures_dir, "Figure_22b.pdf"), 
       smoking_coffee_consumption_direction_plot, 
       width = 10, height = 7)


#' @title Subsection 6.3 - Effect of BMI on the risk of PD
#' @author Ioan Gabriel Bucur


source('R/write_configuration_file.R')
source('R/PolyChord_interface.R')

library(bayesplot)
library(ggplot2)
library(ggthemes)


# Deriving the sufficient statistics --------------------------------------

J <- nrow(genetic_associations)
dat <- genetic_associations[order(genetic_associations$gamma_hat, decreasing = T)[1:J],]

MAF <- dat$theta
gamma_hat <- dat$gamma_hat
gamma_se <- dat$gamma_se
Gamma_hat <- dat$Gamma_hat
Gamma_se <- dat$Gamma_se
beta_hat <- 0 # meta-analysis (Wang)

n <- 2
gamma_N <- 339224
Gamma_N <- 13708 + 95282
beta_N <- 430854 # meta-analysis (Wang)

N <- min(beta_N, gamma_N, Gamma_N)

SSS <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, beta_N, MAF, n)
# write.table(SSS, file = 'data/BMI_PD_SSc.txt', row.names = FALSE, col.names = FALSE)

# for the reverse model, reverse X and Y (last two variables)
RSS <- SSS[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]
# write.table(RSS, file = 'data/BMI_PD_SSr.txt', row.names = FALSE, col.names = FALSE)

# standard deviation of binomial genetic variants
SiG <- sqrt(n * MAF * (1 - MAF))
# write.table(SiG, file = 'data/BMI_PD_SiG.txt', row.names = FALSE, col.names = FALSE)

# Writing the PolyChord configuration file --------------------------------


write_configuration_file(
  config_file_name = 'ini/smmr_6_3.ini',
  SS_filename = 'data/BMI_PD_SSc.txt',
  sigmaG_filename = 'data/BMI_PD_SiG.txt',
  num_instruments = 77,
  num_observations = 108990,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 500)
)

system('bin/polychord_MR ini/smmr_6_3.ini')

beta <- quick_derive_beta(read.table('chains/smmr_6_2_fig_17_equal_weights.txt'), spike = 10000)
beta <- as.matrix(beta)
colnames(beta) <- 'beta'

BW_FG_beta_plot <- mcmc_areas(beta, pars = "beta") +
  theme_tufte() + coord_cartesian(xlim = c(-0.4, 0.2)) + theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.155) + geom_vline(xintercept = c(-0.223, -0.088), linetype = "dashed") + xlab(expression(beta)) + ylab("Posterior density") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -4.3, vjust = -3.1, parse = T, size = 10)

plot(BW_FG_beta_plot)
ggsave('Manuscript-figure17.pdf', BW_FG_beta_plot, width = 10, height = 5)


# Second hyperparameter setting (lambda = 1e-4) ---------------------------

write_configuration_file(
  config_file_name = 'ini/smmr_6_2_fig_18.ini',
  SS_filename = 'data/BW_FG_SSc.txt',
  sigmaG_filename = 'data/BW_FG_SiG.txt',
  num_instruments = 7,
  num_observations = 855,
  spike_precision = 10000,
  model_type = 1,
  PolyChord_control = list(num_live_points = 500)
)

system('bin/polychord_MR ini/smmr_6_2_fig_18.ini')

beta <- quick_derive_beta(read.table('chains/smmr_6_2_fig_18_equal_weights.txt'), spike = 1)
beta <- as.matrix(beta)
colnames(beta) <- 'beta'

BW_FG_beta_plot <- mcmc_areas(beta, pars = "beta") +
  theme_tufte() + coord_cartesian(xlim = c(-0.4, 0.2)) + theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.155) + geom_vline(xintercept = c(-0.223, -0.088), linetype = "dashed") + xlab(expression(beta)) + ylab("Posterior density") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -4.3, vjust = -3.1, parse = T, size = 10)

plot(BW_FG_beta_plot)
ggsave('Manuscript-figure18.pdf', BW_FG_beta_plot, width = 10, height = 5)

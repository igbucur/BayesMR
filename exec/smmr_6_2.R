## Section 6 - Real-world applications
## 6.1 - Effect of birth weight on adult fasting glucose

source('R/write_configuration_file.R')
source('R/PolyChord_interface.R')

library(bayesplot)
library(ggplot2)
library(ggthemes)


# First hyperparameter setting (lambda = 1e-2) ----------------------------

write_configuration_file(
  config_file_name = 'ini/smmr_6_2_fig_17.ini',
  SS_filename = 'data/BW_FG_SSc.txt',
  sigmaG_filename = 'data/BW_FG_SiG.txt',
  num_instruments = 7,
  num_observations = 855,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 500)
)

system('bin/polychord_MR ini/smmr_6_2_fig_17.ini')

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

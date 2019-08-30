#' @title Subsection 6.4 - Does coffee consumption influence smoking?
#' @author Ioan Gabriel Bucur
#' @description Experiment for the purpose of investigating the relationship between
#' coffee consumption and smoking. It is unclear whether this association is causal,
#' and if so which is the cause and which is the effect.
#' @references 


source('R/write_configuration_file.R')
source('R/PolyChord_interface.R')
source('R/derive_sufficient_statistics.R')


# Deriving the sufficient statistics --------------------------------------
n <- 2
J <- 8


EAF <- c(0.61, 0.89, 0.63, 0.29, 0.28, 0.81, 0.33, 0.45)
sigma_G <- sqrt(n * EAF * (1 - EAF))

# Ware et al. 2017 (https://onlinelibrary.wiley.com/doi/pdf/10.1111/add.13888), Table 3 (Page 1849) + Table S3
# Coffee consumption (coffee cups) and smoking heaviness (cigarettes per day) in UK Biobank

gamma_hat <- c(0.03, 0.03, 0.05, 0.06, 0.05, 0.03, 0.09, 0.03)
gamma_se <- c(0.01, 0.02, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01)
gamma_N <- 30062

Gamma_hat <- c(-0.214, 0.203, 0.016, -0.006, 0.015, 0.207, 0.076, -0.144)
Gamma_se <- c(0.132, 0.202, 0.128, 0.140, 0.137, 0.167, 0.135, 0.125)

Gamma_N <- obs_N <- 8072
beta_hat <- 0.45

# varIV <- (Gamma_se / gamma_hat)^2 + (Gamma_hat * gamma_se / gamma_hat / gamma_hat)^2

N <- min(gamma_N, Gamma_N, obs_N)



SSS <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, obs_N, EAF)
RSS <- SSS[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]



# Writing the PolyChord configuration file --------------------------------


write_configuration_file(
  config_file_name = 'ini/smmr_6_4_dir.ini',
  SS_filename = 'data/caffeine_nicotine_SSc.txt',
  sigmaG_filename = 'data/caffeine_nicotine_SiG.txt',
  num_instruments = J,
  num_observations = N,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

# WARNING: The posterior computation takes a lot of time (days) and requires a lot of memory.
system('bin/polychord_MR ini/smmr_6_4_dir.ini')


write_configuration_file(
  config_file_name = 'ini/smmr_6_4_rev.ini',
  SS_filename = 'data/caffeine_nicotine_SSr.txt',
  sigmaG_filename = 'data/caffeine_nicotine_SiG.txt',
  num_instruments = J,
  num_observations = N,
  spike_precision = 10000,
  PolyChord_control = list(num_live_points = 1000)
)

system('bin/polychord_MR ini/smmr_6_4_rev.ini')

# Derive log-evidence and produce plot ------------------------------------

lev_dir <- read_polychord_logevidence('chains/smmr_6_4_dir.stats')
lev_rev <- read_polychord_logevidence('chains/smmr_6_4_rev.stats')
# lev_dir <- -0.422878161448718E+005
# lev_dir_se <- 0.140144123185408E+000
# lev_rev <- -0.422876287259804E+005
# lev_rev_se <- 0.133930400352899E+000

prob_model <- c(1, exp(lev_rev$lev - lev_dir$lev)) / (1 + exp(lev_rev$lev - lev_dir$lev))

beta_dir <- data.frame(quick_derive_beta(read.table('chains/smmr_6_4_dir_equal_weights.txt'), spike = 10000))
colnames(beta_dir) <- 'beta'

dir_p <- ggplot(beta_dir, aes(x = beta)) +
  geom_density(aes_q(y = bquote(..density.. * .(prob_model[1]))), fill = 'gray') +
  ggthemes::theme_tufte(base_size = 30) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  annotate("text", x = -0.25, y = 1.5, label = sprintf("%.2f%%", round(prob_model[2] * 100, 2)), size = 8) +
  xlim(c(-1, 2)) +
  xlab("Causal effect of coffee consumption on smoking") +
  ylab("Posterior density")

ggsave("Manuscript-figure22a.pdf", dir_p, width = 10, height = 7)
# setEPS()
# postscript("../tex/MR Paper/fig/causal_effect_coffee_smoking.eps", fonts = "serif")
# dir_p
# dev.off()

beta_rev <- data.frame(quick_derive_beta(read.table('chains/smmr_6_4_rev_equal_weights.txt'), spike = 10000))
colnames(beta_rev) <- 'beta'

rev_p <- ggplot(beta_rev, aes(x = beta)) +
  geom_density(aes_q(y = bquote(..density.. * .(prob_model[2]))), fill = 'gray') +
  ggthemes::theme_tufte(base_size = 30) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  annotate("text", x = -0.03, y = 60, label = sprintf("%.2f%%", round(prob_model[1] * 100, 2)), size = 8) +
  xlim(c(-0.2, 0.2)) +
  xlab("Causal effect of smoking on coffee consumption") +
  ylab("Posterior density")

ggsave("Manuscript-figure22b.pdf", rev_p, width = 10, height = 7)

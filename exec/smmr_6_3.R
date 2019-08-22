#' @title Subsection 6.3 - Effect of BMI on the risk of PD
#' @author Ioan Gabriel Bucur


source('R/write_configuration_file.R')
source('R/PolyChord_interface.R')

library(bayesplot)
library(ggplot2)
library(ggthemes)


# Deriving the sufficient statistics --------------------------------------

genetic_associations <- read.csv('data/SNP_bmi_parkinson.csv')
J <- nrow(genetic_associations)
dat <- genetic_associations[order(genetic_associations$gamma_hat, decreasing = T)[1:J],]

EAF <- dat$theta
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

# WARNING: The posterior computation takes a lot of time (days) and requires a lot of memory.
system('bin/polychord_MR ini/smmr_6_3.ini')

# As a faster alternative, one can sample from the posterior using Stan,
# since for this example we are not interested in computing the marginal likelihood.
# library(rstan)
# Pks_BMI_fit <- stan("src/BayesMR.stan", data = list(
#   N = N,
#   J = J,
#   n = n,
#   v_spike = 0.0001,
#   v_slab = 1,
#   SS = SSS,
#   mu_X = 0,
#   mu_Y = 0,
#   eaf = EAF
# ), iter = 50000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))
# beta <- as(Pks_BMI_fit, "matrix")[, 'beta', drop = FALSE]


beta <- quick_derive_beta(read.table('chains/smmr_6_3_equal_weights.txt'), spike = 10000)
beta <- as.matrix(beta)
colnames(beta) <- 'beta'

Pks_BMI_beta_plot <- mcmc_areas(beta, pars = "beta") +
  theme_tufte() + theme(text = element_text(size = 30), axis.title.x = element_text(margin = margin(t = 10)), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.1946717) + geom_vline(xintercept = c(-0.3682688, -0.02107458), linetype = "dashed") + xlim(c(-0.5, 0.1)) +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -4.6, vjust = -5, parse = T, size = 8) + ylab("Posterior density") +
  xlab(expression(paste("Effect of BMI on the risk of PD (", beta, ")")))

ggsave('Manuscript-figure19.pdf', Pks_BMI_beta_plot, width = 10, height = 5)


# Reproduce outlier plot --------------------------------------------------

outliers <- order(Gamma_hat)[1:2] # The outliers are the ones with the smallest association with the outcome
colours <- rep('black', J); colours[outliers] <- 'red'
shapes <- rep(1, J); shapes[outliers] <- 2

Pks_BMI_outliers_plot <- ggplot(data.frame(Y = Gamma_hat, X = gamma_hat), aes(x = X, y = Y)) + 
  geom_point(colour = colours, shape = shapes, size = 3) + 
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col = 'red', linetype = 2) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, data = data.frame(Y = Gamma_hat, X = gamma_hat)[-outliers,], col = 'black') +
  theme_tufte() + theme(text = element_text(size = 25), axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 10))) +
  xlab("Genetic associations with the exposure") +
  ylab("Genetic associations with the outcome")

ggsave('Manuscript-figure20.pdf', Pks_BMI_outliers_plot, fonts = "serif", width = 10, height = 7)

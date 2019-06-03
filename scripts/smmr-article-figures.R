#' @author Ioan Gabriel Bucur
#' @description Script to create figures for SMMR Mendelian Randomization article

library(ggplot2)
library(ggthemes)
library(bayesplot)

source("utils/PolyChord_interface.R")

# List of figures ---------------------------------------------------------

# CHECK Figure 14 (Posterior near-LCD)
# CHECK Figure 15 (Estimation Robustness)
# CHECK Figure 16 (Hyperparameter Sensitivity)
# Figure 17 (Causal Effect BW FG)
# Figure 18 (Pleiotropy BW FG)
# Figure 19 (Causal Effect BMI PD)
# Figure 20 (Outliers BMI PD)
# Figure 21 (Pleiotropy BMI PD)
# Figure 22 (Causal Effect Smoking Caffeine)


# Simulation Figures ------------------------------------------------------


model_averaging_plot <- function(root_dir, root_rev, J, dir = TRUE, xlim = c(-0.25, 2)) {
  
  model_comparison <- derive_polychord_Bayes_factor(
    paste0(root_dir, '.stats'), paste0(root_rev, '.stats')
  )
  
  if (dir) {
    PC_dir <- read.table(paste0(root_dir, '_equal_weights.txt'))
    beta_dir <- quick_read_beta(PC_dir, spike = 100)
  } else {
    PC_rev <- read.table(paste0(root_rev, '_equal_weights.txt'))
    beta_rev <- quick_read_beta(PC_rev, spike = 100)
  }
  
  if (dir) {
    fig <- ggplot(as.data.frame(beta_dir), aes(x = beta_dir)) +
      geom_density(aes_q(y = bquote(..density.. * .(model_comparison$prob_dir))), fill = 'gray') +
      theme_tufte(base_size = 30) + xlim(xlim) +
      geom_vline(xintercept = 0, lty = 'dashed') +
      xlab(expression(beta[X %->% Y])) +
      ylab("Posterior density")
  } else {
    fig <- ggplot(as.data.frame(beta_rev), aes(x = beta_rev)) +
      geom_density(aes_q(y = bquote(..density.. * .(model_comparison$prob_rev))), fill = 'gray') +
      theme_tufte(base_size = 30) + xlim(xlim) +
      geom_vline(xintercept = 0, lty = 'dashed') +
      annotate("text", x = annot_x, y = annot_y, label = sprintf("%.2f%%", round(model_comparison$prob_dir * 100, 2))) +
      xlab(expression(beta[Y %->% X])) +
      ylab("Posterior density")
  }
  
  # Do annotation separately
  
  if (dir) {
    fig <- fig + annotate("text", x = xlim[2] / 12, y = 0.9 * max(ggplot_build(fig)$data[[1]]$y), 
                          label = sprintf("%.2f%%", round(model_comparison$prob_rev * 100, 2)), size = 8)
  } else {
    fig <- fig + annotate("text", x = xlim[2] / 12, y = 0.9 * max(ggplot_build(fig)$data[[1]]$y), 
                          label = sprintf("%.2f%%", round(model_comparison$prob_dir * 100, 2)), size = 8)
  }
  
  fig
}


# Near-LCD Example --------------------------------------------------------


nearLCD_causal_effect <- model_averaging_plot(
  'chains/weak_ply_weak_conf/direct/spkrob_PC_dir_sp2_2000_weak',
  'chains/weak_ply_weak_conf/reverse/spkrob_PC_rev_sp2_2000_weak',
  1, xlim = c(-0.25, 2.25)
)

ggsave('../tex/MR Paper/fig/Manuscript-figure14.pdf', nearLCD_causal_effect, width = 10, height = 5)


## Plots associated with BiMR table

# biMR_plots <- list()
# 
# for (se in seq(-0.5, 0.5, 0.1)) {
#   
#   print(paste("se = ", se))
#   
#   biMR_plots[[as.character(se)]] <- model_averaging_plot(
#     sprintf('chains/bidirectional/direct/biMR_PC_dir_sp2_2000_bidirectional_alpha=%+1.2f', se),
#     sprintf('chains/bidirectional/reverse/biMR_PC_rev_sp2_2000_bidirectional_alpha=%+1.2f', se),
#     2, xlim = c(-0.25, 2.25)
#   )
#   ggsave(file = sprintf("~/surfdrive/tex/MR Paper/fig/biMR_PC_dir_sp2_2000_alpha=%+1.2f.pdf", se), 
#          device = "pdf", width = 10, height = 5, dpi = 600)
# }

# Estimation Robustness ---------------------------------------------------

estrob_data <- readRDS('results/iv_validity_robustness.rds')
names(estrob_data) <- c('100%', '80%', '60%', '40%', '20%', '0%')

robust_pleiotropy_posterior <- mcmc_areas(estrob_data) +
  theme_tufte() + coord_cartesian(xlim = c(-0.5, 1.5)) +
  theme(text = element_text(size = 30), axis.title.x = element_text(margin = margin(t = 10))) +
  geom_vline(xintercept = 1) +
  xlab(expression(beta)) +
  ylab("Percentage of valid instruments")

ggsave('../tex/MR Paper/fig/Manuscript-figure15.pdf', robust_pleiotropy_posterior, width = 10, height = 7)



# Sensitivity to Hyperparameters ------------------------------------------

senhyp_data <- readRDS('results/weak_conf_weak_ply_polychord_spike_robustness.rds')

sensitivity_hyperparameters <- ggplot(senhyp_data, aes(x = spike, y = beta)) +
  geom_violin(fill = 'gray', colour = 'darkblue') +
  ylab(expression(beta))+
  xlab(expression(lambda)) +
  coord_flip(ylim = c(-0.5, 1.5)) +
  theme_tufte(base_size = 30) +
  theme(aspect.ratio = 2) +
  theme(axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = sapply(paste0("10^", -as.integer(levels(senhyp_data$spike))), 
                                   function(t) parse(text = t), 
                                   USE.NAMES = FALSE)) +
  scale_linetype_manual(name = "", values = c("WPP bounds" = "dashed")) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  geom_hline(yintercept = 1, lty = 2)

ggsave('../tex/MR Paper/fig/Manuscript-figure16.pdf', sensitivity_hyperparameters, height = 10, width = 6)


# Low birth weight vs Adult fasting gluocse -------------------------------

BW_FG_fit <- readRDS('results/BW_FG_daly_sp4_fixed.rds')
BW_FG_mat <- as.matrix(BW_FG_fit)

BW_FG_beta_plot <- mcmc_areas(BW_FG_mat, pars = "beta") +
  theme_tufte() + coord_cartesian(xlim = c(-0.4, 0.2)) + theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.155) + geom_vline(xintercept = c(-0.223, -0.088), linetype = "dashed") + xlab(expression(beta)) + ylab("Posterior density") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -4.3, vjust = -3.1, parse = T, size = 10)

ggsave('../tex/MR Paper/fig/Manuscript-figure17.pdf', BW_FG_beta_plot, width = 10, height = 5)

BW_FG_fit_test <- readRDS('results/BW_FG_daly_sp4_adapt.rds')
BW_FG_mat_test <- as.matrix(BW_FG_fit_test)
# 
BW_FG_beta_plot <- mcmc_areas(BW_FG_mat_test, pars = "beta") +
  theme_tufte() + coord_cartesian(xlim = c(-0.4, 0.2)) + theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.155) + geom_vline(xintercept = c(-0.223, -0.088), linetype = "dashed") + xlab(expression(beta)) + ylab("Posterior density") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -4.3, vjust = -3.1, parse = T, size = 10)

ggsave('../tex/MR Paper/fig/Manuscript-figure18.pdf', BW_FG_beta_plot, width = 10, height = 5)

# BW_FG_alpha <- mcmc_areas(BW_FG_mat, regex_pars = "^alpha") +
#   theme_tufte() + theme(text = element_text(size = 30)) +
#   #coord_cartesian(xlim = c(-0.08, 0.04)) +
#   scale_y_discrete(labels = parse(text = paste0("alpha[", 1:7, "]")))



# BW_FG_fit_sp4 <- readRDS('results/BW_FG_daly_sp4_fixed.rds')
# BW_FG_mat_sp4 <- as.matrix(BW_FG_fit_sp4)

# BW_FG_beta_plot <- mcmc_areas(cbind(beta = BW_FG_mat[, 'beta'], beta_test = BW_FG_mat_test[, 'beta'])) +
#   theme_tufte() + coord_cartesian(xlim = c(-0.4, 0.2)) + theme(text = element_text(size = 30)) +
#   geom_vline(xintercept = -0.155) + geom_vline(xintercept = c(-0.223, -0.088), linetype = "dashed") +
#   scale_y_discrete(labels = c(1e-2, 1e-4)) + ylab(expression(lambda)) + xlab(expression(beta)) +
#   annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -3.1, vjust = -3.1, parse = T, size = 10)



# Parkinson's versus BMI --------------------------------------------------

load('results/Pks_BMI_sp4.RData')
genetic_associations <- read.csv('data/SNP_bmi_parkinson.csv')

J <- nrow(genetic_associations)
MAF <- genetic_associations$theta
gamma_hat <- genetic_associations$gamma_hat
gamma_se <- genetic_associations$gamma_se
Gamma_hat <- genetic_associations$Gamma_hat
Gamma_se <- genetic_associations$Gamma_se

plot_Pks_BMI_beta <- mcmc_areas(Pks_BMI_x, pars = "beta") +
  theme_tufte() + theme(text = element_text(size = 30), axis.title.x = element_text(margin = margin(t = 10)), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = -0.1946717) + geom_vline(xintercept = c(-0.3682688, -0.02107458), linetype = "dashed") +
  annotate(geom="text", x=-Inf, y=-Inf, label = "hat(beta)^{IVW}", hjust = -5.2, vjust = -10, parse = T, size = 8) +
  ggplot2::scale_y_discrete(labels = "") + ylab("Posterior density") +
  xlab(expression(paste("Effect of BMI on the risk of PD (", beta, ")")))

ggsave('../tex/MR Paper/fig/Manuscript-figure19.pdf', plot_Pks_BMI_beta, width = 8, height = 8)

plot_Pks_BMI_outliers <- ggplot(data.frame(Y = Gamma_hat, X = gamma_hat), aes(x = X, y = Y)) + 
  geom_point(colour = c('red', 'red', rep('black', J - 2)), shape = c(2, 2, rep(1, J - 2)), size = 3) + 
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col = 'red', linetype = 2) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, data = data.frame(Y = Gamma_hat, X = gamma_hat)[-c(1, 2),], col = 'black') +
  theme_tufte() + theme(text = element_text(size = 25), axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 10))) +
  xlab("Genetic associations with the exposure") +
  ylab("Genetic associations with the outcome")

ggsave('../tex/MR Paper/fig/Manuscript-figure20.pdf', plot_Pks_BMI_outliers, fonts = "serif", width = 10, height = 7)

plot_Pks_BMI_alpha <- mcmc_areas(Pks_BMI_x, pars = c("alpha[1]", "alpha[2]")) +
  theme_tufte() + theme(text = element_text(size = 30), axis.title.x = element_text(margin = margin(t = 30))) +
  ggplot2::scale_y_discrete(labels = c("rs17001654", "rs13107325")) +
  xlab(expression(paste("Pleiotropic effect on risk of PD (", alpha, ")")))

ggsave('../tex/MR Paper/fig/Manuscript-figure21.pdf', plot_Pks_BMI_alpha, width = 10, height = 8)


# Coffee Consumption versus Heaviness of Smoking

# evidence when spike = 1e-4, obtained with MN IS 2000 / 5000 live points?
# lev_dir <- -42280.595391969051
# lev_rev <- -42280.297690317493
# lev_dir_se <- 7.5945860427106035E-002
# lev_rev_se <- 7.5551551285966120E-002

# evidence when spike = 1e-4, obtained with PolyChord 2000 live points
lev_dir <- -0.422878161448718E+005
lev_dir_se <- 0.140144123185408E+000
lev_rev <- -0.422876287259804E+005
lev_rev_se <- 0.133930400352899E+000

prob_model <- c(1, exp(lev_rev - lev_dir)) / (1 + exp(lev_rev - lev_dir))

load("results/caffeine_nicotine_sp4_direct.RData")
dir_samples <- as.data.frame(dfit)

load('results/caffeine_nicotine_sp4_reverse.RData')
rev_samples <- as.data.frame(rfit)

dir_p <- ggplot(dir_samples, aes(x = beta)) +
  geom_density(aes_q(y = bquote(..density.. * .(prob_model[1]))), fill = 'gray') +
  theme_tufte(base_size = 30) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  annotate("text", x = -0.25, y = 1.5, label = sprintf("%.2f%%", round(prob_model[2] * 100, 2)), size = 8) +
  xlim(c(-1, 2)) +
  xlab("Causal effect of coffee consumption on smoking") +
  ylab("Posterior density")

ggsave("../tex/MR Paper/fig/Manuscript-figure22a.pdf", dir_p, width = 10, height = 7)
# setEPS()
# postscript("../tex/MR Paper/fig/causal_effect_coffee_smoking.eps", fonts = "serif")
# dir_p
# dev.off()

rev_p <- ggplot(rev_samples, aes(x = beta)) +
  geom_density(aes_q(y = bquote(..density.. * .(prob_model[2]))), fill = 'gray') +
  theme_tufte(base_size = 30) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  annotate("text", x = -0.03, y = 60, label = sprintf("%.2f%%", round(prob_model[1] * 100, 2)), size = 8) +
  xlim(c(-0.2, 0.2)) +
  xlab("Causal effect of smoking on coffee consumption") +
  ylab("Posterior density")

ggsave("../tex/MR Paper/fig/Manuscript-figure22b.pdf", rev_p, width = 10, height = 7)
# setEPS()
# postscript("../tex/MR Paper/fig/causal_effect_smoking_coffee.eps", fonts = "serif")
# rev_p
# dev.off()

# MN_data <- read.table('~/ownCloud/CHiLL/dev/MultiNest_v3.11/chains/sp4_5000_reverse-.txt')
# sigma_X <- qnorm(MN_data[, 22] / 2 + 0.5, sd = 10)
# sigma_Y <- qnorm(MN_data[, 23] / 2 + 0.5, sd = 10)
# # beta_slab <- qnorm(MN_data[, 21]) * sigma_Y / sigma_X
# # beta_spike <- qnorm(MN_data[, 21], sd = 0.01) * sigma_Y / sigma_X
# beta <- sapply(MN_data[, 19], quantile_spike_and_slab_2, spike = 10000) * sigma_Y / sigma_X
# # beta <- sbeta * sigma_Y / sigma_X
# plot(density(beta, weights = MN_data[,1]))
# plot(density(MN_data[, 4], weights = MN_data[,1]))
# 
# DN_data <- read.table('~/dev/DNest4/code/Examples/BayesianMR/posterior_sample.txt')
# DN_w <- read.table('~/dev/DNest4/code/Examples/BayesianMR/weights.txt')
# sigma_X <- qnorm(DN_data[, 22] / 2 + 0.5, sd = 10)
# sigma_Y <- qnorm(DN_data[, 23] / 2 + 0.5, sd = 10)
# sbeta <- sapply(DN_data[, 19], quantile_spike_and_slab_2)
# beta <- sbeta * sigma_Y / sigma_X
# plot(density(beta))





# Produce radial plot
ivw_radial_result <- ivw_radial(format_radial(gamma_hat, Gamma_hat, gamma_se, Gamma_se, RSID = 1:77), weights = 1, alpha = 0.005)
egger_radial_result <- egger_radial(format_radial(gamma_hat, Gamma_hat, gamma_se, Gamma_se, RSID = 1:77), weights = 1, alpha = 0.005)
Pks_BMI_radial_plot <- plot_radial(c(ivw_radial_result, egger_radial_result), F, F, F)
Pks_BMI_radial_plot + theme_tufte() + theme(legend.title = element_blank(), legend.position = "bottom", text = element_text(size = 20))

# Caffeine versus nicotine consumption ------------------------------------

load('results/caffeine_nicotine_sp2_direct.RData')



source('utils/PolyChord_interface.R')
source('MR/run-PolyChord-bidirectional.R')

# Generate data for weak pleiotropy, weak confounding scenario
set.seed(1918)
samples <- 1e4
MAF <- 0.3
n <- 2
G <- rbinom(samples, n, MAF)
U <- rnorm(samples)
X <- G + 0.1 * U + rnorm(samples)
Y <- 0.1 * G + X + 0.1 * U + rnorm(samples)


Z <- cbind(1, G, X, Y)
sigma_G <- sqrt(MAF * (1 - MAF) * n)
SSS <- t(Z) %*% Z / samples
RSS <- SSS[c(1, 2, 4, 3), c(1, 2, 4, 3)]

file_wpwc_SSc <- 'dev/PolyChord/data/weak_ply_weak_conf_SSc.txt'
file_wpwc_SSr <- 'dev/PolyChord/data/weak_ply_weak_conf_SSr.txt'
file_wpwc_SiG <- 'dev/PolyChord/data/weak_ply_weak_conf_SiG.txt'

write.table(SSS, file = file_wpwc_SSc, row.names = FALSE, col.names = FALSE)
write.table(RSS, file = file_wpwc_SSr, row.names = FALSE, col.names = FALSE)
write.table(sigma_G, file = file_wpwc_SiG, row.names = FALSE, col.names = FALSE)

# PolyChord config
n_live <- 2000

for (i in 0:6) {
  spike_precision <- 10^i
  write_polychord_config("dev/PolyChord/ini/config_dir.ini", paste0(sprintf("spkrob_PC_dir_sp%d_", i), n_live, "_weak"), n_live, 9)
  write_model_config("dev/PolyChord/ini/model_dir.ini", 1, 'dev/PolyChord/data/weak_ply_weak_conf_SSc.txt', 'dev/PolyChord/data/weak_ply_weak_conf_SiG.txt', spike_precision = spike_precision)
  
  system('dev/PolyChord/bin/polychord_MR_dir dev/PolyChord/ini/config_dir.ini')
  
  write_polychord_config("dev/PolyChord/ini/config_rev.ini", paste0(sprintf("spkrob_PC_rev_sp%d_", i), n_live, "_weak"), n_live, 9)
  write_model_config("dev/PolyChord/ini/model_rev.ini", 1, 'dev/PolyChord/data/weak_ply_weak_conf_SSr.txt', 'dev/PolyChord/data/weak_ply_weak_conf_SiG.txt', spike_precision = spike_precision)
  
  system('dev/PolyChord/bin/polychord_MR_rev dev/PolyChord/ini/config_rev.ini')
}


Bf_list <- list()

for (i in 0:6) {
  Bf_list[[as.character(i)]] <- derive_polychord_Bayes_factor(
    sprintf("chains/weak_ply_weak_conf/direct/spkrob_PC_dir_sp%d_2000_weak.stats", i),
    sprintf("chains/weak_ply_weak_conf/reverse/spkrob_PC_rev_sp%d_2000_weak.stats", i)
  )
}

PC_beta <- list()

for (i in 0:6) {
  PC_data <- read.table(sprintf("chains/weak_ply_weak_conf/direct/spkrob_PC_dir_sp%d_2000_weak_equal_weights.txt", i))
  PC_beta[[as.character(i)]] <- quick_read_beta(PC_data, spike = 10^i)        
}

df <- data.frame(spike = integer(0), samples = numeric(0))

for (i in 0:6) {
  df <- rbind(df, cbind(rep(i, length(PC_beta[[as.character(i)]])), PC_beta[[as.character(i)]]))
}

names(df) <- c("spike", "beta")
df$spike <- as.factor(df$spike)


p <- ggplot(df, aes(x = spike, y = beta)) +
  geom_violin(fill = 'gray', colour = 'darkblue') +
  ylab(expression(beta))+
  xlab(expression(lambda)) +
  coord_flip(ylim = c(-0.5, 1.5)) +
  theme_tufte(base_size = 30) +
  theme(aspect.ratio = 2) +
  theme(axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = sapply(paste0("10^", -as.integer(levels(df$spike))), 
                                   function(t) parse(text = t), 
                                   USE.NAMES = FALSE)) +
  scale_linetype_manual(name = "", values = c("WPP bounds" = "dashed")) +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  geom_hline(yintercept = 1, lty = 2)


ggsave('../tex/MR Paper/fig/spike_robustness_weak_ply_weak_cnf.pdf', height = 10, width = 6, device = 'pdf')

#' @title Association between Coffee Consumption and Heaviness of Smoking
#' @author Ioan Gabriel Bucur
#' @description Experiment for the purpose of evaluating the effect of coffee consumption on heaviness of smoking


library(rstan)

n <- 2
J <- 8

# note: actually EAF
MAF <- c(0.61, 0.89, 0.63, 0.29, 0.28, 0.81, 0.33, 0.45)
sigma_G <- sqrt(n * MAF * (1 - MAF))

gamma_hat <- c(0.03, 0.03, 0.05, 0.06, 0.05, 0.03, 0.09, 0.03)
gamma_se <- c(0.01, 0.02, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01)
gamma_N <- 30062

# Coffee cups vs Smoking CPD (Ware)
# Gamma_hat <- c(-0.0821, 0.0297, -0.0679, -0.0881, 0.1387, 0.0465, -0.2125, -0.0642)
# Gamma_se <- c(0.0826, 0.1355, 0.0906, 0.0973, 0.1377, 0.1045, 0.0964, 0.0817)
# Gamma_N <- 38181

# Coffee cups vs Smoking cotinine (Ware)
# Gamma_hat <- c(-0.0102, -0.0139, 0.0087, 0.0114, -0.0319, 0.0163, -0.0647, 0.0002)
# Gamma_se <- c(0.0207, 0.0323, 0.0208, 0.0284, 0.0258, 0.0262, 0.0284, 0.0204)

# Coffee Cups and Cigarettes per Day (Ware, Biobank)
Gamma_hat <- c(-0.214, 0.203, 0.016, -0.006, 0.015, 0.207, 0.076, -0.144)
Gamma_se <- c(0.132, 0.202, 0.128, 0.140, 0.137, 0.167, 0.135, 0.125)


Gamma_N <- obs_N <- 8072
beta_hat <- 0.45

# varIV <- (Gamma_se / gamma_hat)^2 + (Gamma_hat * gamma_se / gamma_hat / gamma_hat)^2

N <- min(gamma_N, Gamma_N, obs_N)


source('utils/derive_sufficient_statistics.R')
SSS <- MR_regression_coefficients_to_moments(J, gamma_hat, gamma_se, gamma_N, Gamma_hat, Gamma_se, Gamma_N, beta_hat, obs_N, MAF)
RSS <- SSS[c(1:(J+1), J+3, J+2), c(1:(J+1), J+3, J+2)]

# direct fit
dfit <- stan("src/multiple_IV_norm.stan", data = list(
  N = N,
  J = J,
  n = n,
  v_spike = 1e-3,
  v_slab = 1,
  SS = SSS,
  mu_X = 0,
  mu_Y = 0,
  theta = MAF
), iter = 20000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))
# library(bridgesampling)
dbridge <- bridge_sampler(dfit)
## Spike 1e0
# lp__ -42277.72
# MN: -42321.582188396562      +/-  0.37205284787146692
# IS: -42296.993348395532      +/-   9.8371417887838458E-002
## Spike 1e1
# lp__ -42262.51 +/- 4.37 
# IS: -42281.897856812080      +/-  0.13113817136204867
#     -42289.670303806030      +/-   9.0452986796142124E-002
#     -42289.611084609416      +/-   9.0413971469693613E-002
# .5e -42275.870566560530      +/-  0.17026475123204049
# .5e -42285.394142000347      +/-   8.6094567567147826E-002
## Spike 1e2
# lp__ -42250.12 +/- 4.42      -42301.97
# MN: -42292.252890658769      +/-  0.28347041831663211 
#     -42292.550048936064      +/-  0.14220060156103029
# IS: -42275.780510045479      +/-  0.11827897940242603
#     -42275.994371238674      +/-  0.11851344336448780
#     -42280.757908932013      +/-   7.9987530640087831E-002
# DN: -42298.0205879
#     -42297.6202552
## Spike 1e3
# lp__: -42226.88 +/- 4.44     -42287.19 / -42286.88 (warp3)
# IS: -42273.093799477349      +/-  0.11107876194414266
#     -42274.648767502535      +/-   7.1598839129873082E-002
## Spike 1e4
# lp__ -42222.01 +/- 4.43      -42289.67 / -42289.65 (warp3)
# MN: -42283.436357173268      +/-  0.24872222229410695
# IS: -42281.106155324887      +/-  0.11989214294335614
#     -42281.984569114502      +/-  0.11959979329736224
#     -42280.595391969051      +/-   7.5945860427106035E-002
# PC: -42287.8161448718        +/-  0.140144123185408
## Spike 1e5
# IS: -42298.426473478889      +/-   8.1863572899010847E-002
#     -42297.349478852208      +/-  0.18450819354499293
#     -42299.482505691623      +/-  0.13084136230146670
## Spike 1e6
# lp__ -42209.95 ?
# MN: -42299.188386432077      +/-  0.25467805549551836 


# reverse fit
rfit <- stan("src/multiple_IV_norm.stan", data = list(
  N = N,
  J = J,
  n = n,
  v_spike = 1e-3,
  v_slab = 1,
  SS = RSS,
  mu_X = 0,
  mu_Y = 0,
  theta = MAF
), iter = 20000, chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 20))
rbridge <- bridge_sampler(rfit, rmodel)
## Spike 1e0
# lp__ -42277.48 +/- 4.2
# MN: -42321.373614608114      +/-  0.37146115422229542
# IS: -42297.009078725001      +/-   9.8393007925572984E-002
## Spike 1e1
# lp__ -42262.12 +/- 4.44 
# IS: -42282.011654565344      +/-  0.13121305068584013
#     -42289.556975146588      +/-   9.0418269430633816E-002
#     -42289.633903526177      +/-   9.0500191140408706E-002
# .5e -42275.924063672610      +/-  0.17052006962248431      
# .5e -42285.593624312161      +/-   8.6278867800823736E-002
## Spike 1e2
# lp__ -42249.07 +/- 4.42     -42301.96
# MN: -42292.842566637126      +/-  0.28547534945312025
#     -42293.213197031320      +/-  0.14326395904540459
# IS: -42275.513790679506      +/-  0.11783267102194368
#     -42275.386361796816      +/-  0.11771331614801801
#     -42280.380982547518      +/-   7.9537367192866748E-002
# DN: -42297.8988656
#     -42297.9440012
## Spike 1e3
# lp__ -42222.67 +/- 3.97      -42286.22 / -42286.47 (warp3)
# IS: -42271.918944777404      +/-  0.10926899023178728
#     -42273.619446607561      +/-   7.0818927832490650E-002
## Spike 1e4
# lp__ -42220.00 +/- 4.01      -42289.86 / -42289.86 (warp3)
# MN: -42282.864599770153      +/-  0.24718860481839564
# IS: -42283.448047484941      +/-  0.12191057745948983
#     -42283.830629205549      +/-  0.12462876311542274
#     -42280.297690317493      +/-   7.5551551285966120E-002
# PC: -42287.6287259804        +/-  0.133930400352899
## Spike 1e5
# IS: -42297.255064963007      +/-   9.8541517004737303E-002
#     -42300.165538389876      +/-  0.18748534059085273
#     -42299.147087871577      +/-  0.13193618603433768
## Spike 1e6
# lp__ -42219.06
# MN: -42302.730934534171      +/-  0.27129453112863766



# Bridge Sampling ---------------------------------------------------------

library(rstan)
library(bridgesampling)
# create empty models for bridge sampling

dmodel <- stan("src/multiple_IV_norm.stan", data = list(
  N = N,
  J = J,
  n = n,
  v_spike = 1e-4,
  v_slab = 1,
  SS = SSS,
  mu_X = 0,
  mu_Y = 0,
  theta = MAF
), chains = 0)
dbridge <- bridge_sampler(dfit, dmodel, method = 'warp3')

rmodel <- stan("src/multiple_IV_norm.stan", data = list(
  N = N,
  J = J,
  n = n,
  v_spike = 1e-4,
  v_slab = 1,
  SS = RSS,
  mu_X = 0,
  mu_Y = 0,
  theta = MAF
), chains = 0)
rbridge <- bridge_sampler(rfit, rmodel, method = 'warp3')

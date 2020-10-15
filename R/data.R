#' A dataset containing results for reproducing Figure 14 in Section 4 of the SMMR article.
#'
#' @format A three-level list containing posterior PolyChord samples in the
#' expected and reverse directions for the near LCD example as well as the
#' log-evidence comparison between the two directions.
"near_LCD_example"

#' A dataset containing results for reproducing Figure 15 in Section 5 of the SMMR article.
#'
#' @format A data frame with six columns containing posterior PolyChord beta
#' (causal effect) samples for six different scenarios in terms of IV validity
#' (100%, 80%, 60%, 40%, 20%, 0%) to show the estimation robustness.
"estimation_robustness"

#' A dataset containing results for reproducing Figure 16 in Section 5 of the SMMR article.
#'
#' @format A data frame containing two columns for showing the sensitivity of
#' the estimation wrt hyperparameters, in this case the spike precision.
#' \describe{
#'   \item{spike}{Exponent of spike component precision (10^exponent).}
#'   \item{beta}{Posterior PolyChord beta (causal effect) samples.}
#' }
"sensitivity_to_hyperparameters"

#' A dataset containing results for reproducing Table 1 in Section 5 of the SMMR article.
#'
#' @format A data frame containin the information of Table 1 which shows how
#' robust the model averaging performed by BayesMR to pleiotropy.
#' \describe{
#'   \item{delta}{The size and direction of the pleiotropic effects.}
#'   \item{expected_gamma_hat}{Genetic association with exposure in expected direction.}
#'   \item{expected_gamma_std_err}{Standard error for genetic association with exposure in expected direction.}
#'   \item{expected_Gamma_hat}{Genetic association with outcome in expected direction.}
#'   \item{expected_Gamma_std_err}{Standard error for genetic association with outcome in expected direction.}
#'   \item{reverse_gamma_hat}{Genetic association with exposure in reverse direction.}
#'   \item{reverse_gamma_std_err}{Standard error for genetic association with exposure in reverse direction.}
#'   \item{reverse_Gamma_hat}{Genetic association with outcome in reverse direction.}
#'   \item{reverse_Gamma_std_err}{Standard error for genetic association with outcome in reverse direction.}
#'   \item{probability_expected}{BayesMR probability of expected as opposed to reverse direction.}
#' }
"model_averaging_robustness"


#' A dataset containing results for reproducing Table 2 in Section 5 of the SMMR article.
#'
#' @format A three level list which contains PolyChord posterior samples and
#' estimated log-evidences in both the expected and reverse direction for seven
#' different values of the nonlinearity parameter A.
"sensitivity_to_nonlinearity"

#' A dataset containing results for reproducing Table 3 in Section 5 of the SMMR article.
#'
#' @format A three level list which contains PolyChord posterior samples and
#' estimated log-evidences in both the expected and reverse direction for seven
#' different values of the nonnormality parameter nu.
"sensitivity_to_nonnormality"

#' A dataset containing results for reproducing Figure 17 in Section 6 of the SMMR article.
#'
#' @format Matrix containing posterior samples of the causal effect from exposure
#' to outcome for the birth weight versus adult fasting glucose example. The
#' inference is performed using the default BayesMR settings.
"birth_weight_fasting_glucose_default"

#' A dataset containing results for reproducing Figure 18 in Section 6 of the SMMR article.
#'
#' @format Matrix containing posterior samples of the causal effect from exposure
#' to outcome for the birth weight versus adult fasting glucose example. The
#' inference is performed by incorporating the strict IV assumptions into the
#' BayesMR model (set model == 1 in configuration file).
"birth_weight_fasting_glucose_IV"

#' A dataset containing results for reproducing Figures 19-21 in Section 6 of the SMMR article.
#' 
#' @format Matrix containing posterior samples of the causal effect from exposure
#' to outcome for the risk of Parkinson's disease versus BMI example.
"Parkinson_BMI_posterior"

#' A dataset containing results for reproducing Figure 22a in Section 6 of the SMMR article.
#' 
#' @format Data frame containing posterior samples for the estimated causal effect 
#' of coffee consumption on the heaviness of smoking.
"coffee_consumption_smoking_direction_posterior"

#' A dataset containing results for reproducing Figure 22b in Section 6 of the SMMR article.
#' 
#' @format Data frame containing posterior samples for the estimated causal effect 
#' of smoking on coffee consumption.
"smoking_coffee_consumption_direction_posterior"

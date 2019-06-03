#' Write PolyChord configuration file for BayesMR problem instances
#'
#' @param config_file_name 
#' @param num_parameters 
#' @param num_live_points 
#' @param num_repeats 
#' @param do_clustering 
#' @param grade_frac 
#' @param precision_criterion 
#' @param calc_weighted_posterior 
#' @param calc_equalw_posterior 
#' @param cluster_posterior 
#' @param boost_posterior 
#' @param output_base_dir 
#' @param output_file_root 
#' @param write_resume_file 
#' @param resume_previous_run 
#' @param write_live_points_file 
#' @param write_dead_points_file 
#' @param feedback_level 
#' @param write_param_names_file 
#'
#' @return
#' @export
#'
#' @examples
write_PolyChord_configuration_file <- function(
  config_file_name,
  num_parameters,
  num_live_points = 100,
  num_repeats = num_parameters * 2,
  do_clustering = TRUE,
  grade_frac = 1,
  precision_criterion = 0.001,
  calc_weighted_posterior = TRUE,
  calc_equalw_posterior = TRUE,
  cluster_posterior = FALSE,
  boost_posterior = 5.0,
  output_base_dir = "chains",
  output_file_root = "MR",
  write_resume_file = FALSE,
  resume_previous_run = FALSE,
  write_live_points_file = FALSE,
  write_dead_points_file = FALSE,
  feedback_level = 1,
  write_param_names_file = TRUE
) {
  
  tryCatch({
    fileConn <- file(config_file_name)
    
    # First section: algorithm settings
    algorithm_settings <- paste(
      "[algorithm settings]",
      paste("nlive =", num_live_points),
      paste("num_repeats =", num_repeats),
      paste("do_clustering =", do_clustering),
      paste("grade_frac =", grade_frac),
      paste("precision_criterion =", precision_criterion),
      sep = "\n"
    )
    
    # Second section: posterior settings
    posterior_settings <- paste(
      "[posterior settings]",
      paste("posteriors =", calc_weighted_posterior),
      paste("equals =", calc_equalw_posterior),
      paste("cluster_posteriors =", cluster_posterior),
      paste("boost_posterior =", boost_posterior),
      sep = "\n"
    )
    
    # Third section: output settings
    output_settings <- paste(
      "[output settings]",
      paste("base_dir =", output_base_dir),
      paste("file_root =", output_file_root),
      paste("write_resume =", write_resume_file),
      paste("read_resume =", resume_previous_run),
      paste("write_live =", write_live_points_file),
      paste("write_dead =", write_dead_points_file),
      paste("feedback =", feedback_level),
      paste("write_paramnames =", write_param_names_file),
      sep = "\n"
    )
    
    # Fourth section: prior settings
    prior_settings <- paste(c("[prior settings]", sapply(1:num_parameters, function(i) {
      sprintf("P : p%02d | \\theta_{%02d} | 1 | uniform | 1 | 0.0 1.0", i, i)
    })), collapse = "\n")
    
    # Last section: derived parameter settings
    derived_parameter_settings <- paste(
      "[derived parameter settings]",
      sep = "\n"
    )
    
    writeLines(c(algorithm_settings, "",
                 posterior_settings, "",
                 output_settings, "",
                 prior_settings, "",
                 derived_parameter_settings), fileConn)
    
  }, warning = function(w) {
    message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
    
  }, error = function(e) {
    message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
    
  }, finally = {
    close(fileConn)
  })
  
}
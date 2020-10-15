#' Function for automatically writing BayesMR configuration files for a given problem.
#'
#' @param config_filename String indicating the name of the configuration file.
#' @param SS_filename String indicating the file where the sufficient statistics are saved.
#' @param sigma_G_filename String indicating the file where the instrument frequency statistics are saved.
#' @param num_instruments Integer number of instruments for the given problem.
#' @param num_observations Integer number of observations in the provided sample.
#' @param slab_precision Numeric precision of prior slab component.
#' @param spike_precision Numeric precision of prior spike component.
#' @param model_type Integer indicating the type of prior considered (0 = spike-and-slab, 1 = Gaussian).
#' @param PolyChord_control List of PolyChord control parameters (see details section).
#' 
#' @details PolyChord control parameters
#'
#' @export
#'
#' @examples
#' write_BayesMR_configuration_file(
#'   config_filename = "ini/BayesMR.ini",
#'   SS_filename = "inst/extdata/BayesMR_SS.txt",
#'   sigma_G_filename = "inst/extdata/BayesMR_sigma_G.txt",
#'   num_instruments = 1, num_observations = 1000
#' )
write_BayesMR_configuration_file <- function(
  config_filename,
  SS_filename,
  sigma_G_filename,
  num_instruments,
  num_observations,
  slab_precision = 1,
  spike_precision = 100,
  model_type = 0,
  PolyChord_control = NULL
) {
  
  default_control_parameters <- list(
    num_live_points = 200,
    num_parameters = num_instruments * 2 + 7,
    num_repeats = 3 * (num_instruments * 2 + 7),
    do_clustering = TRUE,
    grade_frac = 1,
    precision_criterion = 0.001,
    calc_weighted_posterior = TRUE,
    calc_equalw_posterior = TRUE,
    cluster_posterior = FALSE,
    boost_posterior = 5.0,
    output_base_dir = "chains",
    output_file_root = tools::file_path_sans_ext(basename(config_filename)),
    write_resume_file = FALSE,
    resume_previous_run = FALSE,
    write_live_points_file = FALSE,
    write_dead_points_file = FALSE,
    feedback_level = 1,
    write_param_names_file = TRUE
  )
  
  # update default control parameters with given input
  for (param in names(PolyChord_control)) {
    default_control_parameters[[param]] <- PolyChord_control[[param]]
  }
  
  tryCatch({
    fileConn <- file(config_filename)
    
    # First section: algorithm settings
    algorithm_settings <- paste(
      "[algorithm settings]",
      paste("nlive =", as.integer(default_control_parameters$num_live_points)),
      paste("num_repeats =", as.integer(default_control_parameters$num_repeats)),
      paste("do_clustering =", default_control_parameters$do_clustering),
      paste("grade_frac =", default_control_parameters$grade_frac),
      paste("precision_criterion =", default_control_parameters$precision_criterion),
      sep = "\n"
    )
    
    # Second section: posterior settings
    posterior_settings <- paste(
      "[posterior settings]",
      paste("posteriors =", default_control_parameters$calc_weighted_posterior),
      paste("equals =", default_control_parameters$calc_equalw_posterior),
      paste("cluster_posteriors =", default_control_parameters$cluster_posterior),
      paste("boost_posterior =", default_control_parameters$boost_posterior),
      sep = "\n"
    )
    
    # Third section: output settings
    output_settings <- paste(
      "[output settings]",
      paste("base_dir =", default_control_parameters$output_base_dir),
      paste("file_root =", default_control_parameters$output_file_root),
      paste("write_resume =", default_control_parameters$write_resume_file),
      paste("read_resume =", default_control_parameters$resume_previous_run),
      paste("write_live =", default_control_parameters$write_live_points_file),
      paste("write_dead =", default_control_parameters$write_dead_points_file),
      paste("feedback =", default_control_parameters$feedback_level),
      paste("write_paramnames =", default_control_parameters$write_param_names_file),
      sep = "\n"
    )
    
    # Fourth section: prior settings
    prior_settings <- paste(c("[prior settings]", sapply(1:default_control_parameters$num_parameters, function(i) {
      sprintf("P : p%02d | \\theta_{%02d} | 1 | uniform | 1 | 0.0 1.0", i, i)
    })), collapse = "\n")
    
    # Fifth section: derived parameter settings
    derived_parameter_settings <- paste(
      "[derived parameter settings]",
      sep = "\n"
    )
    
    # Last section: model settings
    
    model_settings <- paste(
      "[model settings]",
      paste("SS_filename =", SS_filename),
      paste("sigma_G_filename =", sigma_G_filename),
      paste("slab_precision =", slab_precision),
      paste("spike_precision =", spike_precision),
      paste("instruments =", as.integer(num_instruments)),
      paste("observations =", as.integer(num_observations)),
      paste("model =", model_type),
      sep = "\n"
    )
    
    writeLines(c(algorithm_settings, "",
                 posterior_settings, "",
                 output_settings, "",
                 prior_settings, "",
                 derived_parameter_settings, "",
                 model_settings), fileConn)
    
  }, warning = function(w) {
    message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
    
  }, error = function(e) {
    message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
    
  }, finally = {
    close(fileConn)
  })
  

}
write_model_configuration_file <- function(
  config_file_name,
  SS_filename,
  sigmaG_filename,
  num_instruments,
  num_observations,
  slab_precision = 1,
  spike_precision = 100,
  model_type = 0
) {
  
  tryCatch({
    fileConn <- file(config_file_name)
    
    model_settings <- paste(
      "[model settings]",
      paste("SS_filename =", SS_filename),
      paste("sigmaG_filename =", sigmaG_filename),
      paste("slab_precision =", slab_precision),
      paste("spike_precision =", spike_precision),
      paste("instruments =", num_instruments),
      paste("observations =", num_observations),
      paste("model =", model_type),
      sep = "\n"
    )
    
    write(model_settings, fileConn)
    
  }, warning = function(w) {
    message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
    
  }, error = function(e) {
    message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
    
  }, finally = {
    close(fileConn)
  })
  

}
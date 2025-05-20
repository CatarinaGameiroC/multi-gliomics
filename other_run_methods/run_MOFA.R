library("MOFA2")
library("reticulate")


run_MOFA <- function(X, num.factors = NULL){

  time_taken <- system.time({
    
    seeds <- seq(1, 10)
    mofa_models <- list()
      
    for (seed in 1: length(seeds)){
  
      mofa <- create_mofa(lapply(X, t))
      
      data_options <- get_default_data_options(mofa)
      model_options <- get_default_model_options(mofa) 
      model_options$spikeslab_weights <- TRUE # (default)
      model_options$ard_weights <- TRUE # (default)

  
      training_options <- get_default_training_options(mofa)
      training_options$seed <- seed
      training_options$convergence_mode <- "medium" # 0.00005% deltaELBO change
    
      if (is.null(num.factors) == FALSE){
        model_options$num_factors <- num.factors
      }
      else {
        training_options$drop_factor_threshold <- 0.02 # factors explaining less than 2% will be dropped
      }
    
      mofa <- prepare_mofa(
        mofa, 
        data_options = data_options, 
        model_options =  model_options, 
        training_options = training_options)
    
      use_condaenv("mofa_env", required = TRUE)
      mofa <- run_mofa(mofa, use_basilisk = FALSE)
      
      mofa_models <- append(mofa_models, mofa)
    }
    
    mofa <- select_model(mofa_models, plot = FALSE)
  
    loadings <- get_weights(mofa)
    factors <- get_factors(mofa)[[1]]
    result <- list(Z = factors, W = loadings, model_output = mofa)
  })
  
  result$time_taken <- time_taken["elapsed"]
  return (result)
}
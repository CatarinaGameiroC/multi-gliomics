library("iClusterPlus")
# there are 3 functions of iCluster in this library ('iCluster', 'iCluster2' and 'iClusterBayes')
# we are going to use the first for now, but the second allows for other types of regularization and clustering,
# while the last, allows for other types of data (and tuning of number of compoenents/clusters)


run_iCluster <- function(X, num.factors = NULL){ 
# INPUTS:
  # 'X': list of omics (each one assumed to be tagged, with rows and columns names, about the same patients), where rows are samples and columns features
  # 'num.factors': number of factors (and the number of clusters is num.factors + 1)
# OUTPUTS:
  # list of two elements: the factor matrix and a list of the loadings matrices for each omics

    
  time_taken <- system.time({
    
    methods <- rep(list(c("lasso")), length(X))
    lasso.penalties_2test <- c(0.01, 0.1, 0.2, 0.3)
    best_pod <- Inf 
    best_model <- NULL
    
    
    if (is.null(num.factors)){
      
      num.factors_2test <- c(2:15)
      
      for (f in num.factors_2test) {
        for (l in lasso.penalties_2test) {
          fit <- iCluster2(x = X, K = f + 1, lambda = lapply(1: length(X), function(x) c(l)), method = methods) 
          pod_value <- compute.pod(fit)
          if (pod_value < best_pod) {
            best_pod <- pod_value
            best_model <- fit
            num.factors <- f
            best_lasso <- l
          # QUESTION: what if the difference is so small that is not siginficant to take into account more number of factors?
          }
        }
      }
    }
    else {
      
      for (l in lasso.penalties_2test){
        fit <- iCluster2(x = X, K =  num.factors + 1, lambda = lapply(1: length(X), function(x) c(l)), method = methods) 
        pod_value <- compute.pod(fit)
        if (pod_value < best_pod) {
          best_pod <- pod_value
          best_model <- fit
          best_lasso <- l
        }
        
      }
    }
    
    icluster <- best_model

    factors <- as.matrix(t(icluster$meanZ))
    rownames(factors) <- rownames(X[[1]])
    colnames(factors) <- 1: num.factors

    loadings <- list()
    for (m in 1: length(X)){
        loadings[[m]] <- as.matrix(icluster$beta[[m]])
        rownames(loadings[[m]]) <- colnames(X[[m]])
        colnames(loadings[[m]]) <- 1: num.factors
    }
    res <- list(Z = factors, W = loadings, clusters = as.matrix(icluster$clusters))
    
  })
  print(paste0("iCluster completed in ", time_taken["elapsed"], " seconds."))
  print(paste0("Number of factors: ", num.factors, " and lambda penalty used: ", best_lasso))
  
  return (res)
}
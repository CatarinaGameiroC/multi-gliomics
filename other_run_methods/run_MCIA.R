library("omicade4")


run_MCIA <- function(X, num.factors = NULL){
# INPUTS:
  # 'X': list of omics (each one assumed to be tagged, with rows and columns names, about the same patients), where rows are samples and columns features
  # 'num.factors': number of factors
# OUTPUTS:
  # list of two elements: the factor matrix and a list of the loadings matrices for each omics

  time_taken <- system.time({
      
    if (is.null(num.factors)){
      mcia <- mcia(lapply(X, t), cia.scan = TRUE)
      num.factors <- mcia$mcoa$nf
    }
    else {
      mcia <- mcia(lapply(X, t), cia.nf = num.factors)
    }
   
    factors <- as.matrix(mcia$mcoa$SynVar)
    rownames(factors) <- rownames(X[[1]])
    colnames(factors) <- 1: num.factors

    loadings <- list()
    ng <- sapply(X, function(x) dim(x)[2])
    ngsum <- cumsum(c(0, ng))
    for (m in 1: length(X)){
      loadings[[m]] <- as.matrix(mcia$mcoa$axis[(ngsum[m] + 1): ngsum[m + 1], ]) 
      rownames(loadings[[m]]) <- colnames(X[[m]])
      colnames(loadings[[m]]) <- 1: num.factors
    }
    res <- list(Z = factors, W = loadings)
  })
  print(paste0("MCIA completed in ", time_taken["elapsed"], " seconds."))
  
  return (res)
}

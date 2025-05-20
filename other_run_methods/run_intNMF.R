library("IntNMF")



run_intNMF <- function(X, num.factors = NULL){
# INPUTS:
  # 'X': list of omics (each one assumed to be tagged, with rows and columns names, about the same patients), where rows are samples and columns features
  # 'num.factors': number of factors, which is in fact the number of clusters assumed 
# OUTPUTS:
  # list of two elements: the factor matrix and a list of the loadings matrices for each omics
    
  time_taken <- system.time({
  
    X_pos <- list()
    for (m in 1: length(X)){
      if (min(X[[m]]) < 0){
          X_pos[[m]] <- X[[m]] + abs(min(X[[m]]))
      } 
      else {
          X_pos[[m]] <- X[[m]]
      }
      X_pos[[m]] <- X_pos[[m]]/ max(X_pos[[m]])
    }
    
    # i think this way has a way of taking into account everything: number of features and missing values.
    denominator <- c()
    for (view in X_pos){
      denominator <- c(denominator, mean(norm(view, type = "F")))
    }
    numerator <- max(denominator)
    wt <- numerator * 1/denominator
      
    if (is.null(num.factors)){

      intnmf <-  nmf.opt.k(dat = X_pos, k.range = 2:15,  make.plot = FALSE, result = TRUE, progress = FALSE,
                            n.runs = 30, n.fold = 5, wt = wt, st.count = 10)
      num.factors <- as.numeric(which.max(apply(intnmf, 1, mean)))

      # NOTE: 'n.runs', 'n.fold', 'wt' and 'st.count' are the default values. 
      # EXPLORE: 'wt' is the importance of each omics: maybe we could adjust this vector to the number of features, 
    # the variation, the missing data of each omics and assign weights according to that (or biological or experimental knowledge)
    # or even use the expression suggested in the paper of the method (that i do not understand).
      # QUESTION: scaling the data before giving it to the model, the weights will have importance right? The omics may have
    # a bif difference on the number of features, e.g..

    }

    intnmf <- nmf.mnnals(dat = X_pos, k = num.factors,
                          maxiter = 200, st.count = 20, n.ini = 30, wt = wt) 

    # NOTE: 'maxiter', 'st.count', 'n.ini' and 'wt' are the default values. 
    # EXPLORE: 'wt' is the importance of each omics: maybe we could adjust this vector to the number of features (in proportion), 
  # the variation represented, the missing data of each omics and assign weights according to that (or biological or experimental knowledge)
  # or even use the expression suggestes in the paper of the method.

    factors <- as.matrix(intnmf$W)
    colnames(factors) <- 1: num.factors
    rownames(factors) <- rownames(X[[1]])

    loadings <- list()
    for (m in 1: length(X)){
      loadings[[m]] <- t(as.matrix(intnmf$H[[m]])) 
      rownames(loadings[[m]]) <- colnames(X[[m]])
      colnames(loadings[[m]]) <- 1: num.factors
    }
  
    res <- list(Z = factors, W = loadings, clusters = as.matrix(intnmf$clusters))
  })
    
  print(paste0("intNMF completed in ", time_taken["elapsed"], " seconds."))
  return (res)
}

# GENERAL NOTES:
    # the method outputs clusters assignments based on: intnmf$clusters == apply(intnmf$Z, 1, which.max) (TRUE in all entries)
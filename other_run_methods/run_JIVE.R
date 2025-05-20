library("r.jive")


run_JIVE <- function(X, num.factors = NULL){
# INPUTS:
  # 'X': list of omics (each one assumed to be tagged, with rows and columns names, about the same patients), where rows are samples and columns features
  # 'num.factors': number of factors
# OUTPUTS:
  # list of two elements: the factor matrix and a list of the loadings matrices for each omics
  
  time_taken <- system.time({
      
    if (is.null(num.factors)){
      jive <- jive(lapply(X, t), showProgress = FALSE, method = "perm", conv = 1e-06, 
                    maxiter = 100, scale = FALSE, center = FALSE, orthIndiv = TRUE, est = TRUE)
    }

    else {
      jive <- jive(lapply(X, t), rankJ = num.factors, rankA = rep(num.factors, length(X)), showProgress = FALSE,
                    method = "given", conv = "default", maxiter = 100, scale = FALSE, center = FALSE, orthIndiv = TRUE, est = TRUE)
    }

    # NOTE: This function centers and scales (dividing by the frobenius norm) the data and replaces the missing values using the 'SVDmiss' function from its library. 
  # 'scale', 'center', 'orthIndiv' and 'est'(compress data using SVD to speed the computational time) are the default values.

    rankJV <- jive$rankJ # ranks of the joint structure of the data (3*features x samples)
    rankIV <- jive$rankA # ranks of the individual structure of the data 

    JV <- numeric(0)
    ng <- 0
    for (m in 1: length(X)){
      JV <- rbind(JV, jive$joint[[m]])
      ng <- c(ng, dim(jive$joint[[m]])[1]) # number of features of each omics
    }
    svd.JV <- svd(JV) # J = U D V^T (=> J^T = V D U^T), the factor matrix will be the truncated (VD) (= S) and the U^T the weight matrices combined of all omics
    S <- svd.JV$v %*% diag(svd.JV$d)
    factors <-  S[ , 1: rankJV]
    rownames(factors) <- rownames(X[[1]])
    colnames(factors) <- 1: rankJV
    
    
    loadings <- list()
    ngsum <- cumsum(ng)
    for (m in 1: length(X)){
      loadings[[m]] <- svd.JV$u[(ngsum[m] + 1) : ngsum[m + 1], 1: rankJV] 
      rownames(loadings[[m]]) <- colnames(X[[m]])
      colnames(loadings[[m]]) <- 1 : rankJV
    }
    res <- list(Z = factors, W = loadings)
  })
  print(paste0("JIVE completed in ", time_taken["elapsed"], " seconds."))

  return (res)
}
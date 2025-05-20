library("RGCCA")


run_RGCCA <- function(X, num.factors = NULL){
  
# INPUTS:
  # 'X': list of omics (each one assumed to be tagged, with rows and columns names, about the same patients), where rows are samples and columns features
  # 'num.factors': number of factors
# OUTPUTS:
  # list of two elements: the factor matrix and a list of the loadings matrices for each omics

  time_taken <- system.time({
    
    # connectivity matrix
    design <- matrix(0.1, ncol = length(X), nrow = length(X))
    diag(design) <- 0
    
    if (is.null(num.factors)){
      
      ncomp_test <- matrix(rep(seq(2,12), length(X)), byrow = FALSE, ncol = length(X))
      num.factors <- rgcca_permutation(X, connection = design, par_type = "ncomp", par_value = ncomp_test, n_perms = 20,
                                        scale = FALSE, method = "rgcca", scheme = "horst", init = "svd", 
                                        verbose = FALSE, NA_method = "na.ignore", quiet = FALSE, n_iter_max = 1000, 
                                        comp_orth = TRUE, scale_block = TRUE, tol = 1e-8)$best_params
    }
    
    # sparsity
    sparsity_matrix <- matrix(rep(c(0.2, 0.5, 0.7), length(X)), byrow = FALSE, ncol = length(X))
    sparsity_opt <- rgcca_permutation(X, connection = design, par_type = "sparsity", par_value = sparsity_matrix, n_perms = 20,
                                        method = "rgcca", scheme = "horst", init = "svd", scale = FALSE,
                                        verbose = FALSE, NA_method = "na.ignore", quiet = FALSE, n_iter_max = 1000, 
                                        comp_orth = TRUE, scale_block = TRUE, tol = 1e-8)$best_params
    
    
    rgcca <- rgcca(X, ncomp = num.factors, connection = design,
                    scheme = "horst", sparsity = sparsity_opt, tau = "optimal", scale = FALSE, 
                    init = "svd", tol = 1e-8, verbose = FALSE,  NA_method = "na.ignore", 
                    quiet = FALSE, n_iter_max = 1000, comp_orth = TRUE, method = "rgcca", scale_block = TRUE)
    
    # centroid: g(x) = |x|
    
    # NOTE: there are three hyper-parameters: 'tau', 'connectivity' and 'sparsity'. the connection matrix is used as the best found in DIABLO. 'sparsity' can be tuned using 'rgcca_permutation', 
  # with 'rgcca_cv' only in a supervised way. 'tau' the same but here is already being optimized.
    
    
    factors <- rgcca$Y 
  
    loadings <- list()
    for (m in 1: length(X)){
      loadings[[m]] <- as.matrix(rgcca$a[[m]])
      rownames(loadings[[m]]) <- colnames(X[[m]])
      colnames(loadings[[m]]) <- 1: num.factors
    }
    res <- list(Z = factors, W = loadings)
  })
  print(paste0("RGCCA completed in ", time_taken["elapsed"], " seconds."))
  
  return (res)
}


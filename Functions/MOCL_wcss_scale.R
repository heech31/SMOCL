
MOCL_wcss_scale <- function(Xf, K, max.iter = 500, df = 3, degree = 2,
                   lambda_seq = seq(0.1, 2, length=50), tol = 1e-6){
  
  ######################
  Xf <- scale(Xf)
  ######################
  n <- nrow(Xf)
  nlambdas <- length(lambda_seq)
  
  clusters_mat <- matrix(0, n, nlambdas)
  lambda_sel_mat <- matrix(0, 5, nlambdas)
  
  colnames(clusters_mat) <- colnames(lambda_sel_mat) <- lambda_seq
  rownames(lambda_sel_mat) <- c("Monotonicity","No 0","wss", "penalty", "nvars")
  
  function_mat_lst <- list()
  coef_mat_lst <- list()
  iter_lst <- list()
  
  lam_mono_check <- 0
  mono_n <- 0
  lam_q <- 1
  
  # If maximum of lambda (max.L) is too small, we generate lambda sequence with 2*max.L value.
  while(lam_mono_check==0){
    if(lam_q>1){ # If monotonicity is not satisfied
      lambda_seq <- lambda_seq + max(lambda_seq)
    }
    
    for(l in 1:length(lambda_seq)){#  l <- 4
      lambda_l <- lambda_seq[l]
      mocl_l <- MOCL_scale(Xf = Xf, K = K, max.iter = max.iter, df = df, degree = degree,
                           lambda = lambda_l, tol = tol)
      
      nvars <- mocl_l$nvars
      penalty <- mocl_l$coop_norm
      wcss <- mocl_l$WCSS
      monotone_flag <- mono_check(mocl_l$ms_coef, mocl_l$groups, ncol(Xf))
      non_null_flag <- zero_check(mocl_l$f_x)
      
      # Save
      clusters_mat[,l] <- mocl_l$cluster
      lambda_sel_mat[,l] <- c(monotone_flag, non_null_flag, wcss, penalty, nvars)
      function_mat_lst[[l]] <- mocl_l$f_x
      iter_lst[[l]] <- mocl_l$iteration
      coef_mat_lst[[l]] <- mocl_l$ms_coef

    }
    
    lam_mono_check = sum(lambda_sel_mat[1,]) # If all lambda values do not guarantee monotonocity -> lam_mono_check=0
    lam_q = lam_q+1
    
  }
  
  # Change the colnames to updated lambda values
  colnames(clusters_mat) <- colnames(lambda_sel_mat) <- lambda_seq
  # Best lambda
  ## lambda satisfying sign-coherence and non-zero
  lambda_good_ind <- which( (lambda_sel_mat[1,]*lambda_sel_mat[2,]) == 1 )
  bics <- (lambda_sel_mat[3,lambda_good_ind] + log(n)*lambda_sel_mat[5,lambda_good_ind])
  best_bic <- min( bics )
  best_lam_ind <- which( bics == best_bic )
    
  if(length(best_lam_ind)!=0){
    best_lam_ind <- min(best_lam_ind)
  }
  ## results for optimal lambda
  lam_wcss <- lambda_seq[ best_lam_ind ]
  cluster_wcss <- clusters_mat[,best_lam_ind]
  fx_wcss <- function_mat_lst[[best_lam_ind]]
  coef_wcss <- coef_mat_lst[[best_lam_ind]]
  
  if(lam_q>2){
    print("Since lambda values are too small to graurantee monotonicity, the values of lambda seqence were enlarged.")
  }
  
  return(
    list(cluster = cluster_wcss, fx = fx_wcss, coef = coef_wcss, lam_wcss = lam_wcss, 
         lambda_sel_mat = lambda_sel_mat,
         clusters_mat = clusters_mat,
         function_mat_lst = function_mat_lst,
         iter_lst = iter_lst)
  )
}




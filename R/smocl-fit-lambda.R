##################################################################################################################
# Monotone Clustering function
##################################################################################################################
MOCL_scale <- function(Xf, K, max.iter = 500, degree = degree, df = 3, lambda = 0.01, tol = 1e-6){
  
  Xf <- (Xf); n <- nrow(Xf); p <- ncol(Xf);
  data_x <- as.matrix(Xf)
  iter <- max.iter

  ## Get spline basis functions
  Z <- monotone.splines(data_x, degree = degree, df = df) 
  if(df >1){
    groups <- rep(1:p, each = df)
  }else if(df==1){
    groups <- (1:p)
  }

  ## Trace objective function
  obj_trace <- matrix(0, max.iter, 4)
  absdiff <- rep(0, max.iter)
  
  ## Initialize
  cluster_matrix <- matrix(0,nrow(data_x),iter+1) # matrix to save cluster rank
  cluster_matrix_order <- matrix(0,nrow(data_x),iter+1) # matrix to save cluster score
  fx_lst <- list()
  coef_lst <- list()
  y_hat_lst <- list()
  y_hat_unscaled_lst <- list()
  
  stop <- 0 # to check the number of iteration
  
  ## Get initials
  initials <- initialize_mocl(data_x, K = K)  
  s <- initials$y
  init_clu_rank <- initials$init_clu_rank
  clu_rank <- init_clu_rank
  ms_coef <- rep(0, ncol(Z))
  

  for(i in 1:iter){# i <- 1

    #################################################################
    # MAM step
    #################################################################
    s_old <- s; ms_coef_old <- ms_coef; clu_rank_old <- clu_rank
    
    ms_lasso <- MS_lasso(Xf = data_x, Zf = Z, Yf = s_old, df = df, lambda = lambda)
    
    if( all(ms_lasso$coef==0) ){
      y_hat <- initials$PC1
      cluster <- initials$init_clu_rank
      break
    }
    
    ms_coef <- ms_lasso$coef
    y_mean <- ms_lasso$b0
    # Yc <- ms_lasso$Yc
    ms_coef <- ms_coef/norm2( y_mean + Z %*% ms_coef )
    
    # prediction
    y_hat = y_mean + Z %*%  ms_coef #+y_mean # estimated clusters score for all observation

    # function estimate
    function_mat <- matrix(0,nrow(data_x),ncol(data_x))
    for(j in 1:ncol(data_x)){
      # coef_j = ms_coef[groups==j]
      coef_j = ms_coef[groups==j]
      z_j = Z[,groups==j, drop=FALSE]
      function_mat[,j] = z_j%*%coef_j
    }
    fx_lst[[i]] = function_mat
    coef_lst[[i]] = ms_lasso$coef
    y_hat_lst[[i]] = y_hat

    #################################################################
    
    
    
    #################################################################
    # UCF step (K-means)
    #################################################################
    
    # For initial centers of UCF step
    ini_centers=rep(0, K)
    for(r in 1:K){
      ini_centers[r] = mean( y_hat[clu_rank==r,] )
    }
    
    km_new <- tryCatch( kmeans(y_hat, centers = c(ini_centers), iter.max=500), error = function(e) 1)
    ucf_count <- 0
    
    while( length(km_new)==1 & ucf_count < 50){# When empty cluster occurs
      km_new <- ClusterR::KMeans_rcpp(y_hat, clusters = K, CENTROIDS = matrix(ini_centers, ncol=1), seed = sample(1e4,1) )
      km_new$centers <- km_new$centroids
      ucf_count <- ucf_count + 1
    }
    
    if (ucf_count >= 50 && length(km_new) == 1) {
      stop(
        "The UCF step failed to obtain a valid K-means solution after 50 attempts.",
        call. = FALSE
      )
    }
    
    o_c <- km_new$cluster # cluster assignments
    
    # relabel index
    clu_rank <- cluster_reassign(clusters = km_new$cluster, centers = km_new$centers) 
    
    # Upadte score vector Uc = s
    oc_order <- o_c*0 
    # oc_order_ini <- o_c*0 
    for(ks in 1:K){
      oc_order[ clu_rank == ks ] <- mean(y_hat[clu_rank==ks])
      # oc_order_ini[ clu_rank_old == ks ] <- mean(y_hat[clu_rank_old==ks])
    }

    cluster_matrix[,i] <- clu_rank # save cluster
    cluster_matrix_order[,i] <- oc_order # save cluster score vector
    
    s_mid <- ini_centers[clu_rank_old]
    s_new <- oc_order
    
    obj_trace[i,1] <- gam_obj(s_old, Z, ms_coef, lambda, p, df)$with_pen
    obj_trace[i,2] <- gam_obj(s_mid, Z, ms_coef, lambda, p, df)$with_pen
    obj_trace[i,3] <- gam_obj(s_new, Z, ms_coef, lambda, p, df)$with_pen
    
    
    absdiff[i] <- abs( obj_trace[i,1] - obj_trace[i,3] )
    if( absdiff[i] < tol){break}
    #update cluster score
    s <- s_new
    # # save number of iteration
    # stop=stop+1
  }
  
  # Get WCSS
  WCSS <- 0
  for( kk in 1:K){
    # sum(y_hat^2)
    WCSS_k <- sum ( ( y_hat[clu_rank == kk] - mean(y_hat[clu_rank == kk]) )^2 )
    WCSS <- WCSS + WCSS_k
  }
  
  
  
  ########################### # retrun
  stop <- i
  
  nvars <- sum(colSums(abs(matrix(ms_coef, df, p)))!=0)
  cluster_matrix <- cluster_matrix[,c(1:(stop))]
  cluster_matrix_score <- cluster_matrix_order[,c(1:(stop))]
  cluster <- clu_rank
  
  if(i==1){
    function_mat <- matrix(0,nrow(data_x),ncol(data_x))
    # for(j in 1:ncol(data_x)){
    #   function_mat[,j] <- pca_dt$rotation[j,1]*data_x[,j]
    # }
    cluster <- init_clu_rank
    cluster_matrix <- cluster
    y_hat <- initials$PC1
  }
  
  mono <- mono_check( ms_coef = ms_coef, ms_group = groups, p = ncol(data_x) )
  coop_norm <- norm_coop(ms_coef, df)
  
  if(stop>=max.iter){
    
    print("no convergence")
    
    convergence = 0
    
    return(list(cluster = cluster, cluster_matrix = cluster_matrix, y_hat = y_hat, mono = mono,
                iteration = stop, f_x = function_mat, fx_lst = fx_lst, coef_lst = coef_lst, y_hat_lst = y_hat_lst,
                y_tr_mean = y_mean, Z = Z, WCSS = WCSS, coop_norm = coop_norm,
                ms_coef = ms_coef, groups = groups, 
                nvars = nvars,
                absdiff = absdiff[1:stop], obj_trace = obj_trace[1:stop,], convergence = convergence) )
    
    
  }else{
    
    convergence = 1
    
    return(list(cluster = cluster, cluster_matrix = cluster_matrix, y_hat = y_hat, mono = mono,
                iteration = stop, f_x = function_mat, fx_lst = fx_lst, coef_lst = coef_lst, y_hat_lst = y_hat_lst,
                y_tr_mean = y_mean, Z = Z, WCSS = WCSS, coop_norm = coop_norm,
                ms_coef = ms_coef, groups = groups, 
                nvars = nvars,
                absdiff = absdiff[1:stop], obj_trace = obj_trace[1:stop,], convergence = convergence) )
  }
  
  
}


#  The "scoop" package from "http://www.math-evry.cnrs.fr/logiciels/scoop" is necessary 


###################################################
# Get rho sequences for simulation
###################################################

get_rho_seq <- function(chr_setting, n_rhos){
  
  if( chr_setting==1 ){
    rho_rel_nom <- c(0.2000000, 0.4494949, 0.6333333, 0.7646465, 0.8828283, 0.9878788, 1.0929293, 1.1717172, 1.2505051, 1.3292929, 1.4080808)
    rho_rel_ord <- rep(1, 11)
  }else{
    rho_rel_nom <- rep(1, 11)
    rho_rel_ord <- c(5.0000000, 2.2318182, 1.5727273, 1.3090909, 1.1333333, 1.0015152, 0.9136364, 0.8696970, 0.7818182, 0.7378788, 0.6939394)
  }
  
  return(list(rho_rel_nom = rho_rel_nom, rho_rel_ord = rho_rel_ord))
  
}


###################################################
# Get MS basis functions
###################################################

monotone.splines <- function(Xf,  degree = 2, df = 3){
  
  if( df == 1){
    
    Z <- Xf
    
  }else{
    Z <- c()
    for( jj in 1:ncol(Xf)){# jj <- 1
      tmpZ <- splines2::iSpline(Xf[,jj], degree = degree, df = df, intercept = FALSE, Boundary.knots = c(min(Xf[,jj])-1.1,max(Xf[,jj]) + 1.1) )
      Z <- cbind(Z, tmpZ)
    }
  }
  
  
  Z <- scale(Z, center = TRUE, scale = TRUE)
  
  return(Z)
  
}


###################################################
# Check sign-coherency
###################################################

mono_check<-function(ms_coef,ms_group,p){
  mono_n=0
  for(j in 1:p){
    coef_g <- ms_coef[ms_group==j]
    n_coef_p <- sum(coef_g>0)
    n_coef_m <- sum(coef_g<0)
    mono_n <- mono_n + (n_coef_p*n_coef_m)
  }
  
  mono <- (mono_n==0)
  
  return(mono)
  
}




###################################################
# Check if at least one non-null function
###################################################
zero_check <- function(function_mat){
  
  no_zero <- (sum(abs(function_mat))>0) # logical
  
  return(no_zero)
}


###################################################
# Monotone Spline Lasso
###################################################

MS_lasso <- function(Xf, Zf, Yf, df, lambda){
  
  p <- ncol(Xf)

  groups <-  rep(1:p, each = df)
  
  coopfit_model <- coop.lasso(Zf, Yf, groups, family = "gaussian", lambda = lambda,
                              intercept = TRUE, normalize = FALSE, eps = 1e-8)
  
  b0 <- slot(coopfit_model,"coefficients")[1]
  coef <- slot(coopfit_model,"coefficients")[-1]
  
  return(list(coef = coef, b0 = b0))
  
}






###############################################################################
# Cluster relabeling
###############################################################################

cluster_reassign <- function(clusters, centers){
  
  K <- max(clusters)
  membership <- clusters
  
  order_scores <- order(centers)
  
  for( kk in 1:K){
    membership[clusters == order_scores[kk]] <- kk
  }
  
  return(membership)
  
}




###############################################################
# minmax scaler for matrix
###############################################################

nor_minmax = function(x){
  xx=matrix(0,nrow(x),ncol(x))
  colnames(xx)=colnames(x)
  for(i in 1:ncol(x)){
    xx[,i]=(x[,i] - min(x[,i])) / (max(x[,i]) - min(x[,i]))
  }
  return(xx)
}


###############################################################
# 2-norm
###############################################################
norm2 <- function(x){
  
  sqrt( sum(x*x) )
  
}


###############################################################
# coop-norm
###############################################################
norm_coop <- function(coef, df){
  
  p <- length(coef)/(df)
  
  groups <- rep(1:p, each = df) # group list for splines
  
  penalty <- 0
  
  for( gg in 1:p){
    
    coef_p <- coef_m <- coef[groups==gg]
    
    coef_p[coef_p<0] <- 0  
    
    coef_m[coef_m>0] <- 0
    
    penalty <- penalty + (norm2(coef_p) + norm2(-coef_m))
    
  }
  
  return(penalty)
  
}



###############################################################
# objective function
###############################################################
gam_obj <- function(s, Z, coef, lambda, p, df){
  
  y_hat <-  Z%*%coef # norm2(y_hat) should be 1
  
  penalty <- norm_coop(coef, df)
  
  without_pen <- 0.5*sum( (s - y_hat)^2 )
  
  # sqrt(df) is the adaptive weight that coop.lasso
  with_pen <- without_pen + lambda*penalty * sqrt(df)
  
  return(list(with_pen = with_pen, without_pen = without_pen, penalty = penalty)  )
  
}
############################################################################


############################################################################
# Get initials for MOCL
############################################################################
initialize_mocl <- function(data, K){
  
  xc <- scale(data, scale=FALSE)
  PC1 <- -(xc %*% svd(xc)$v[,1])[,,drop=TRUE]
  PC1 <- PC1 / norm2( PC1 )
  
  km_PC1 <- kmeans(as.matrix(PC1), center = K, nstart = 100, iter.max = 100 )
  
  yhat_init <-  PC1 # initial estimated cluster score (continuous)
  yhat_init <-  yhat_init
  init_clu <- km_PC1$cluster # initial cluster assignments
  y <-  yhat_init*0  # for discrete initial cluster score 
  
  for(ks in 1:K){
    
    y[init_clu==ks]=mean(yhat_init[init_clu==ks])
    
  }

  init_clu_rank <- cluster_reassign(km_PC1$cluster, km_PC1$centers) # rank of initial ordinal clusters
  
  res <- list(PC1 = PC1, y = y,  init_clu_rank = init_clu_rank)
  
}



############################################################################
group_means <- function( xs, ys, group){

  xmeans <- aggregate(xs ~ group, FUN = mean)
  ymeans <- aggregate(ys ~ group, FUN = mean)

  res <- cbind( xmeans[,2], ymeans[,2] )

  return( res )

}


# ############################################################################
# # CH index
# ############################################################################
# CHindex <- function(x, clusters){
#   
#   x <- scale(x, scale=FALSE)
#   n <- length(clusters)
#   K <- max(clusters)
#   U <- matrix(0, n, K)
#   
#   for( ii in 1:n){
#     
#     U[ii,clusters[ii]] <- 1
#     
#   }
#   
#   sizes <- colSums(U)
#   xcc <- (U%*%(crossprod(U, x)/sizes))
#   
#   num <- sum(colSums((xcc)^2))/(K-1)
#   denom <- sum(colSums((x-xcc)^2))/(n-K)
#   
#   return(num/denom)
# }
# ############################################################################
# 
# 
# ############################################################################
# # Get between and within cluster sum of squares
# get_css <- function(X, clusters){
#   K <- max(clusters)
#   z <- matrix(0, nrow(X), K)
#   for( kk in 1:K){
#     z[clusters==kk,kk] <- 1
#   }
#   
#   prj_center <- z%*%solve(crossprod(z),t(z))
#   means <- crossprod(z,X)/colSums(z)
#   
#   bcss <- sum(diag(t(X)%*%prj_center%*%X))
#   wcss <- sum(diag(t(X)%*%(diag(nrow(X)) - prj_center)%*%X))
#   
#   return(list(bcss = bcss, wcss = wcss))
# }
# ############################################################################
# 
# 
# 

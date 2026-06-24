
gendata <- function(rho_ord = 1, rho_nom = 1, cov_setting = 1,
                    ord_size = c(50,50,50,50), 
                    nom_size = c(50,50,50,50), 
                    s_ord = 3, s_nom = 4, p = 20){
  
  # cov_setting == 1 => all identity
  # cov_stting == 2  => random cov (see gen_covs.R)
  
  K_ord <- length(ord_size)
  K_nom <- length(nom_size)
  n <- sum(ord_size)
  pi_k_ord <- ord_size/n
  pi_k_nom <- nom_size/n
  
  clu_ord <- rep(1:K_ord, times = ord_size)
  clu_nom <- rep(1:K_nom, times = nom_size)
  perm <- sample(n, n, replace=FALSE)
  clu_nom <- clu_nom[perm]
  

  mu_mat <- matrix(0, s_ord, K_ord)
  mu_mat[1,] <-  rho_ord * 10/(1+exp(-3*(1:K_ord-2.5)))
  mu_mat[2,] <- -rho_ord * (exp(1:K_ord)/5)
  mu_mat[3,] <-  rho_ord * (9*log(1:K_ord))
  
  # Data generation
  
  ############################
  # Ordinal clusters    ######
  ############################
  if( cov_setting == 1){
    Covs <- rep( list(diag(s_ord)), K)
  }else{
    Covs <- mapply(gen_covs, rep(s_ord,K), SIMPLIFY = FALSE)
  }
  
  x <- c()
  for( kk in 1:K){
    x <- rbind(x, rmvnorm(sum(clu_ord==kk), mean = mu_mat[,kk], sigma= Covs[[kk]]) )
  }
  x_ord <-  x
  colnames(x_ord) <- paste("ord",1:(ncol(x_ord )),sep="_")
  
  pop_bcss_ord <- sum( rowSums( scale(t(mu_mat), scale=FALSE)^2 ) * pi_k_ord )
  
  pop_wcss_ord <- sum( unlist( lapply(Covs, function(x) sum(diag(x)) ) ) * pi_k_ord )
  
  ch_ord <- pop_bcss_ord/pop_wcss_ord
  
  ############################
  # Nominal clusters    ######
  ############################
  if( cov_setting == 1){
    Covs <- rep( list(diag(s_nom)), K)
  }else{
    Covs <- mapply(gen_covs, rep(s_nom,K), SIMPLIFY = FALSE)
  }
  
  mu_no_mat = diag(K)
  mu_no = 10
  # Nominal cluster generation
  x_nom <-  matrix( rnorm(n*ncol(mu_no_mat)), n, ncol(mu_no_mat))

  for( kk in 1:K){# kk <- 1
    x_nom[,kk] <- rho_nom * mu_no * (clu_nom == kk) + x_nom[,kk]
  }
  colnames(x_nom) <- paste("nom",1:(ncol(x_nom )),sep="_")
  
  mu_mat <- rho_nom * mu_no * mu_no_mat
  pop_bcss_nom <- sum( rowSums( scale(t(mu_mat), scale=FALSE)^2 ) * pi_k_nom )
  
  pop_wcss_nom <- sum( unlist( lapply(Covs, function(x) sum(diag(x)) ) ) * pi_k_nom )
  
  ch_nom <- pop_bcss_nom/pop_wcss_nom
  
  ###############################
  # Noise variables
  p_e <- p - s_ord - s_nom
  
  x_noise <- matrix(rnorm(n*p_e,0,1),n,p_e)
  colnames(x_noise) <- paste("n",1:(ncol(x_noise )),sep="_")
  
  dat <- scale( cbind(x_ord, x_nom, x_noise) )
  
  
  return( list( dat = dat, 
                clu_ord = clu_ord, clu_nom = clu_nom,
                ch_nom = ch_nom, ch_ord = ch_ord) )
  
  
}


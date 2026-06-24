gen_covs <- function(s){
  
  evalues <- 10^seq(0,-2,l=s)
  
  tmp_svd <- svd( matrix(rnorm(s^2),s,s) )
  
  R <- tmp_svd$u %*% t(tmp_svd$v)
  
  Cov <- R %*% (t(R) * (evalues))
  
  Cor <- cov2cor( Cov )
  
  return(Cor)
  
}


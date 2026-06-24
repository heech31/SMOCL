stability_K <- function(data1, data2, sample_id1, sample_id2, 
                              sample_id_intersect, lambda_seq){
  
  res1 <- MOCL_wcss_scale(Xf = data1,  K = K,  lambda_seq  = lambda_seq, df = 5)$cluster
  res2 <- MOCL_wcss_scale(Xf = data2,  K = K,  lambda_seq  = lambda_seq, df = 5)$cluster
  

  # clusters1
  cluster_ind_1 <- lapply(res1, function(x) x[ match(sample_id_intersect, sample_id1)])
  cluster_ind_2 <- lapply(res2, function(x) x[ match(sample_id_intersect, sample_id2)])
  
  cors <- mapply(function(c1,c2) cor(c1, c2, method = "kendall"), cluster_ind_1, cluster_ind_2)
  airs <- mapply(function(c1,c2) pdfCluster::adj.rand.index(c1, c2), cluster_ind_1, cluster_ind_2)
  
  return( data.frame(cors = cors, airs = airs) )

}



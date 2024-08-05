library(sparcl)

Sparse_km<-function(data,kk=kk,wbounds_list=NULL,nperms_n=10){
  s_km_per=KMeansSparseCluster.permute(data,K=kk,wbounds=wbounds_list,nperms=nperms_n,silent = TRUE)
  bw=(s_km_per$wbounds==s_km_per$bestw)
  gap_1sd=s_km_per$gaps[bw]-s_km_per$sdgaps[bw] # 1sd gap
  final_w=s_km_per$wbounds[s_km_per$gaps>gap_1sd][1] # s based on 1sd gap
  
  # s_kmeans with wbound based on 1sd
  km.out_best<-KMeansSparseCluster(data,K=kk,wbounds=s_km_per$bestw)
  sparse_m_best=t(as.matrix(((km.out_best[[1]]$ws)^2>0)*1))
  weigth_best=km.out_best[[1]]$ws
  
  km.out_1sd <- KMeansSparseCluster(data,K=kk,wbounds=final_w)
  sparse_m_1sd=t(as.matrix(((km.out_1sd[[1]]$ws)^2>0)*1))
  weigth_1sd=km.out_1sd[[1]]$ws 
  
  return(list(sparse_1sd=sparse_m_1sd,cluster_1sd=km.out_1sd[[1]]$Cs,w_1sd=final_w,weigth_1sd=weigth_1sd,
              sparse_best=sparse_m_best,cluster_best=km.out_best[[1]]$Cs,w_best=s_km_per$bestw,weigth_best=weigth_best))
}

getAnywhere(KMeansSparseCluster.permute)

# Auto MOCL
MOCL_wcss<-function(Xf,k,max.iter=500,num.knotsf=6,
                   lambda_seq=seq.default(from=0.001,to=1,length=50),seed=100,delta=0.01,r_a=0.01){
  
  clusters_mat=matrix(0,nrow(Xf),length(lambda_seq))
  lambda_sel_mat=matrix(0,3,length(lambda_seq))
  colnames(clusters_mat)=colnames(lambda_sel_mat)=lambda_seq
  rownames(lambda_sel_mat)=c("Monotonicity","No 0","wss")
  function_mat_lst = list()
  coef_mat_lst=list()
  iter_lst = list()
  
  lam_mono_check=0
  mono_n=0
  lam_q=1
  
  # If maximum of lambda (max.L) is too small, we generate lambda sequence with 2*max.L value.
  while(lam_mono_check==0){
    if(lam_q>1){ # If monotonicity is not satisfied
      lambda_seq=lambda_seq+max(lambda_seq)
    }
    
    for(l in 1:length(lambda_seq)){
      lambda_l=lambda_seq[l]
      mocl_l=MOCL(Xf=Xf,k=k,max.iter=max.iter,num.knotsf=num.knotsf,lambda=lambda_l,delta=delta)
    
      lam_sel=lambda_selection(ms_Z=mocl_l$Z,ms_yhat=mocl_l$y_hat,ms_Yc=mocl_l$Yc,ms_coef=mocl_l$ms_coef,
                               function_mat=mocl_l$f_x,
                               ms_group=mocl_l$groups,r_a=r_a,ms_cluster=mocl_l$cluster)
      # Save
      clusters_mat[,l]=mocl_l$cluster
      lambda_sel_mat[,l]=c(lam_sel$mono,lam_sel$no_zero,lam_sel$wss_ratio)
      function_mat_lst[[l]]=mocl_l$f_x
      iter_lst[[l]]=mocl_l$iteration
      coef_mat_lst[[l]]=mocl_l$ms_coef

    }
    
    lam_mono_check=sum(lambda_sel_mat[1,]) # If all lambda values do not guarantee monotonocity -> lam_mono_check=0
    lam_q=lam_q+1
  }
  
  # change the colnames to updated lambda values
  colnames(clusters_mat)=colnames(lambda_sel_mat)=lambda_seq
  # Best lambda
  ## lambda satisfying sign-coherence and no 0 condition
  lambda_sel_mat_opt=lambda_sel_mat[,(lambda_sel_mat[1,]*lambda_sel_mat[2,]==1)]
  lambda_seq_opt=lambda_seq[(lambda_sel_mat[1,]*lambda_sel_mat[2,]==1)]
  
  ## optimal lambda corresponding to wss ratio
  lam_wcss=lambda_seq_opt[which.min(lambda_sel_mat_opt[3,])]
  wcss_w=which(lambda_seq==lam_wcss)
  
  
  ## results for optimal lambda
  cluster_wcss=clusters_mat[,wcss_w]
  fx_wcss=function_mat_lst[[wcss_w]]
  coef_wcss=coef_mat_lst[[wcss_w]]
  
  if(lam_q>2){
    print("Since lambda values are too small to graurantee monotonicity, the values of lambda seqence were enlarged.")
  }
  
  return(list(cluster=cluster_wcss,fx=fx_wcss,coef=coef_wcss,lam_wcss=lam_wcss,
              lambda_sel_mat=lambda_sel_mat,clusters_mat=clusters_mat,function_mat_lst=function_mat_lst,iter_lst=iter_lst))
}



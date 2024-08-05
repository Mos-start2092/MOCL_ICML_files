rand_summary<-function(mm_lst){
  for(i in 1:length(mm_lst)){
    if(i==1){
      rand_m=mm_lst[[i]]$rand_mat
    }else{rand_m=rbind(rand_m,mm_lst[[i]]$rand_mat)}
  }
  rand_mean=apply(rand_m,2,mean)#mean)
  rand_2se=2*apply(rand_m,2,sd)/sqrt(length(mm_lst))
  return(list(rand_m=rand_m,rand_mean=rand_mean,rand_2se=rand_2se))
}

rand_no_summary<-function(mm_lst){
  
  for(i in 1:length(mm_lst)){
    mm=mm_lst[[i]]
    
    rand_no_lst=rep(0,nrow(mm$cluster_mat))
    
    for(j in 1:nrow(mm$cluster_mat)){
      rand_no_lst[j]=rand.index(mm$clu_no,mm$cluster_mat[j,])
    }
    
    if(i==1){
      rand_no_mat=t(matrix(c(rand_no_lst)))
    }else{
      rand_no_mat=rbind(rand_no_mat,c(rand_no_lst))
    }   
  }
  colnames(rand_no_mat)=rownames(mm$cluster_mat)
  rand_no_mean=apply(rand_no_mat,2,mean)
  
  return(list(rand_no_mat=rand_no_mat,rand_no_mean=rand_no_mean))
}

adj_rand_summary<-function(mm_lst){
  for(i in 1:length(mm_lst)){
    if(i==1){
      rand_m=mm_lst[[i]]$adj_rand_mat
    }else{rand_m=rbind(rand_m,mm_lst[[i]]$adj_rand_mat)}
  }
  rand_mean=apply(rand_m,2,mean)
  rand_2se=2*apply(rand_m,2,sd)/sqrt(length(mm_lst))
  return(list(adj_rand_m=rand_m,adj_rand_mean=rand_mean,adj_rand_2se=rand_2se))
}

kendall_summary<-function(mm_lst){
  for(i in 1:length(mm_lst)){
    if(i==1){
      summary_m=mm_lst[[i]]$kendall_mat
    }else{summary_m=rbind(summary_m,mm_lst[[i]]$kendall_mat)}
  }
  summary_mean=apply(abs(summary_m),2,mean)
  summary_2se=2*apply(abs(summary_m),2,sd)/sqrt(length(mm_lst))
  return(list(kendall_m=summary_m,kendall_mean=summary_mean,kendall_2se=summary_2se))
}

TP_FP_summary<-function(mm_lst){
  for(i in 1:length(mm_lst)){
    if(i==1){
      TP_m=mm_lst[[i]]$TP_FP_mat[,1]
      FP_m=mm_lst[[i]]$TP_FP_mat[,2]
    }else{
      TP_m=rbind(TP_m,mm_lst[[i]]$TP_FP_mat[,1])
      FP_m=rbind(FP_m,mm_lst[[i]]$TP_FP_mat[,2])
    }
  }
  precision_m=TP_m/(TP_m+FP_m)
  return(list(TP_m=TP_m,FP_m=FP_m,precision_m=precision_m,
              TP_mean=apply(TP_m,2,mean),FP_mean=apply(FP_m,2,mean),precision_mean=apply(precision_m,2,mean)))
}

norm_vec <- function(x) sqrt(sum(x^2))

coef_summary<-function(mm_lst){
  knot_s=nrow(mm_lst[[1]]$coef_mocl)/ncol(mm_lst[[1]]$data)
  
  for(i in 1:length(mm_lst)){
    mm=mm_lst[[i]]
    kendall_mat=data.frame(mm$kendall_mat)
    sign_mat=data.frame(2*(kendall_mat>0)-1)
    
    coef_mocl_i<-as.numeric(sign_mat["MOCL_wss"])*mm$coef_mocl
    coef_admocl_i<-as.numeric(sign_mat["Ad_MOCL"])*mm$coef_admocl
    
    coef_mocl_i=coef_mocl_i/norm_vec(c(coef_mocl_i))
    coef_admocl_i=coef_admocl_i/norm_vec(c(coef_admocl_i))
    
    coef_mocl_i_mean=apply(matrix(coef_mocl_i,8),2,mean)
    coef_admocl_i_mean=apply(matrix(coef_admocl_i,8),2,mean)
    
    if(i==1){
      coef_mocl_m=coef_mocl_i_mean
      coef_admocl_m=coef_admocl_i_mean
    }else{
      coef_mocl_m=rbind(coef_mocl_m,coef_mocl_i_mean)
      coef_admocl_m=rbind(coef_admocl_m,coef_admocl_i_mean)
    }
  }
  return(list(coef_mocl_m=coef_mocl_m,coef_admocl_m=coef_admocl_m))
  
}

coef_mocl_summary<-function(mm_lst){
  knot_s=nrow(mm_lst[[1]]$coef_mocl)/ncol(mm_lst[[1]]$data)
  
  for(i in 1:length(mm_lst)){
    mm=mm_lst[[i]]
    sign=mm$kendall/abs(mm$kendall)
    
    coef_mocl_i<-as.numeric(sign)*mm$coef_mocl
    
    coef_mocl_i=coef_mocl_i/norm_vec(c(coef_mocl_i))
    
    coef_mocl_i_mean=apply(matrix(coef_mocl_i,8),2,mean)

    if(i==1){
      coef_mocl_m=coef_mocl_i_mean
    }else{
      coef_mocl_m=rbind(coef_mocl_m,coef_mocl_i_mean)
    }
  }
  return(list(coef_mocl_m=coef_mocl_m))
}

weight_summary<-function(mm_lst){
  for(i in 1:length(mm_lst)){
    mm=mm_lst[[i]]
    wei_skm_m=mm$weight_skm/sd(c(mm$weight_skm))
    wei_skm_1sd_m=mm$weight_skm_1sd/sd(c(mm$weight_skm_1sd))
    if(i==1){
      weight_skm_m=wei_skm_m
      weight_skm_1sd_m=wei_skm_1sd_m
    }else{
      weight_skm_m=rbind(weight_skm_m,wei_skm_m)
      weight_skm_1sd_m=rbind(weight_skm_1sd_m,wei_skm_1sd_m)
    }
  }
  colnames(weight_skm_m)=colnames(weight_skm_1sd_m)=c(1:ncol(mm$data))
  return(list(weight_skm_m=weight_skm_m,weight_skm_1sd_m=weight_skm_1sd_m))
}

precision_recall<-function(models_l,true_index=c(1,2,3)){
  TP=apply(models_l$sparse_mat[,true_index],1,sum)
  FP=apply(models_l$sparse_mat[,-true_index],1,sum)
  FN=apply((1-models_l$sparse_mat[,true_index])^2,1,sum)
  TN=apply((1-models_l$sparse_mat[,-true_index])^2,1,sum)
  
  precision=TP/(TP+FP)
  recall=TP/(TP+FN)
  f1=2*recall*precision/(recall+precision)
  return(list(precision=precision,recall=recall,f1=f1))
}

pre_re_f1_summary<-function(mm_lst,true_index=c(1,2,3)){
  for(i in 1:length(mm_lst)){
    p_r=precision_recall(mm_lst[[i]],true_index=true_index)
    if(i==1){
      precision_m=p_r$precision
      recall_m=p_r$recall
      f1_m=p_r$f1
    }else{
      precision_m=rbind(precision_m,p_r$precision)
      recall_m=rbind(recall_m,p_r$recall)
      f1_m=rbind(f1_m,p_r$f1)
    }
  }
  precision_m[is.na(precision_m)]<-0
  recall_m[is.na(recall_m)]<-0
  f1_m[is.na(f1_m)]<-0
  
  precision_mean=round(apply(precision_m,2,mean),3)
  precision_2se=round(apply(precision_m,2,sd)/sqrt(i)*2,3)
  
  recall_mean=round(apply(recall_m,2,mean),3)
  recall_2se=round(apply(recall_m,2,sd)/sqrt(i)*2,3)
  
  f1_mean=round(apply(f1_m,2,mean),3)
  f1_2se=round(apply(f1_m,2,sd)/sqrt(i)*2,3)
  
  return(list(precision_mean=precision_mean,precision_2se=precision_2se,
              recall_mean=recall_mean,recall_2se=recall_2se,
              f1_mean=f1_mean,f1_2se=f1_2se,f1_mat=f1_m))
}




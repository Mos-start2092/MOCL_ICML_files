# Functions for lambda selection

###################################################
###################################################
# gaurantee monotonocity
mono_check<-function(ms_coef,ms_group,p){
  mono_n=0
  for(j in 1:p){
    coef_g=ms_coef[ms_group==j]
    n_coef_p=sum(coef_g>0)
    n_coef_m=sum(coef_g<0)
    mono_n=mono_n+(n_coef_p*n_coef_m)
  }
  mono=(mono_n==0)
  return(mono)
}
###################################################
###################################################
# no 0
zero_check<-function(function_mat){
  no_zero=(sum(abs(function_mat))>0)
  return(no_zero)
}

########################################################################################


WSS_TSS<-function(ms_yhat,ms_cluster){
  kk=unique(ms_cluster)
  tss<-sum((ms_yhat-mean(ms_yhat))^2)
  wss=0
  for(k in 1:kk){
    y_k=ms_yhat[ms_cluster==k]
    wss=wss+sum((y_k-mean(y_k))^2)
  }
  
  wss_ratio=wss/tss

  return(wss_ratio=wss_ratio)
}


########################################################################################
########################################################################################
lambda_selection<-function(ms_Z,ms_Yc,ms_coef,function_mat,ms_group,r_a,ms_yhat,ms_cluster){
  p=length(unique(ms_group))
  
  mono=mono_check(ms_coef=ms_coef,ms_group=ms_group,p=p)
  no_zero=zero_check(function_mat=function_mat)
  wss_ratio=WSS_TSS(ms_yhat=ms_yhat,ms_cluster=ms_cluster)
  
  return(list(wss_ratio=wss_ratio,mono=mono,no_zero=no_zero))
}

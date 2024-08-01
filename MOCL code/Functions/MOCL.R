#  The "scoop" package from "http://www.math-evry.cnrs.fr/logiciels/scoop" is necessary 
library(scoop)

# Monotone spline basis generator
## functions for computation of I-Splines of degree 2,cf Ramsay 88
## It's from "mslasso2.r" package

isplineDesign <- function(x,knot.vec){
  I <- numeric()
  deg <- 2
  for(i in 1:(length(knot.vec)-deg)){
    if(x<knot.vec[i]){
      I[i] <- 0
    }
    if (x>knot.vec[i+2]) {
      I[i] <- 1
    }else{
      if(x>=knot.vec[i] && x<=knot.vec[i+1]){
        I[i] <- (x-knot.vec[i])^2 / ((knot.vec[i+1]-knot.vec[i])*(knot.vec[i+2]-knot.vec[i]))
      }
      if(x>=knot.vec[i+1] && x<=knot.vec[i+2]){
        I[i] <- 1-(knot.vec[i+2]-x)^2 / ((knot.vec[i+2]-knot.vec[i])*(knot.vec[i+2]-knot.vec[i+1]))
      }
    }
  }
  return(I)
}

monotone.basis <- function(X,k,spline,xl,xr,intercept){
  n <- length(X)
  dx <- (xr-xl)/(k-1)
  num.knots <- k-2
  
  if(spline=="ispline"){
    deg <- 2
    knot.vec <- seq(xl-(deg-1)*dx,xr+(deg-1)*dx,by=dx) 
    knot.vec[(deg+1):(num.knots+2)] <- quantile(unique(X), seq(0,1,length = (num.knots+2)))[-c(1,(num.knots+2))]
    B <- t(apply(as.matrix(X),1,"isplineDesign",knot.vec=knot.vec))
  } 
  res <- list("B"=B, "knot.vec"=knot.vec)
  return(res)
}

monotone.splines <- function(Xf,num.knots){
  pf <- ncol(Xf)
  Z <- NULL
  for(j in 1:pf){
    x <- Xf[,j]
    spline <- monotone.basis(x, (num.knots+2), "ispline", xl = min(x), xr = max(x), FALSE)
    z1 <- spline$B
    Z=cbind(Z,z1)
  }
  Z <- scale(Z, center = TRUE, scale = FALSE)
  return(Z)
}

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

########################################################################################
MS_lasso<-function(Xf,Zf,Yf,num.knotsf,lambda){
  groups <- as.vector(t(matrix(rep((1):(ncol(Xf)),(num.knotsf+2)),ncol(Xf),(num.knotsf+2)))) # group list for splines
  Yc<-Yf-mean(Yf) # centered cluster score
  b0=mean(Yf)
  lambda_s=lambda
  coopfit_model <- coop.lasso(Zf,Yc, groups, family = "gaussian", intercept = FALSE, lambda = lambda_s)
  coef=slot(coopfit_model,"coefficients")
  
  return(list(coef=coef,b0=b0,Yc=Yc))
}
  
  
  
#######################################################################################


library(ClusterR) # For K-means++

###############################################################################
# reassign cluster function -> to reorder clusters
## we give clusters the order corresponding to y values.
cluster_reassign<-function(y,new_cluster,k){
  original_cluster=y
  n=length(original_cluster)
  c_mat=cbind(original_cluster,new_cluster)
  
  # reassign by mean of original_cluster up to new_cluster
  cc_mat=matrix(0,k,2)
  cc_mat[,1]=c(1:k)
  
  for(i in 1:k){
    cc_mat[i,2]=mean(c_mat[(c_mat[,2]==i),1])
  }
  
  colnames(cc_mat)=c("new_cluster","mean_original_data")
  cc_mat=data.frame(cc_mat)
  cc_mat=cc_mat[order(cc_mat$mean_original_data),]
  
  new_cc=c(1:k)
  cc_mat=cbind(cc_mat,new_cc)
  cc_mat=cc_mat[,c(1,3)]     # remove mean
  
  # merge with new_cluster(non_ordinal_cluster)
  index=c(1:n)
  ordinal_clu_mat=cbind(new_cluster,index)
  ordinal_clu_mat=merge(ordinal_clu_mat,cc_mat,by="new_cluster")
  ordinal_clu_mat=ordinal_clu_mat[order(ordinal_clu_mat$index),]
  
  # new_ordinal_cluster
  new_ordinal_cluster=ordinal_clu_mat$new_cc
  return(new_ordinal_cluster)
}

###############################################################
# minmax scaler for matrix
nor_minmax = function(x){
  xx=matrix(0,nrow(x),ncol(x))
  colnames(xx)=colnames(x)
  for(i in 1:ncol(x)){
    xx[,i]=(x[,i] - min(x[,i])) / (max(x[,i]) - min(x[,i]))
  }
  return(xx)
}

##################################################################################################################
# Monotone Clustering function
## MSLasso_bic.R function is necessary.
MOCL<-function(Xf,k,max.iter=500,num.knotsf=6,lambda=0.1,delta=0.01){
  
  Xf=nor_minmax(Xf);n=nrow(Xf);p=ncol(Xf);data_x=as.matrix(Xf)
  iter=max.iter
  kk=k
  lambda=lambda
  
  #########################################################################
  #########################################################################
  # PC1 -> K-means++ -> initial cluster score
  pca_dt <- prcomp(data_x,center = T,scale. = T) # scaling -> PCA
  PC1=pca_dt$x[,1] # PC1
  km_PC1=KMeans_rcpp(as.matrix(PC1), clusters = kk, num_init = 5, max_iters = 200, 
                     initializer = 'kmeans++')### apply K-means++ to PC1
  
  yhat_init=PC1 # initial estimated cluster score (continuous)
  init_clu=km_PC1$cluster # initial cluster assignments
  y=yhat_init*0  # for discrete initial cluster score 
  for(ks in 1:kk){
    y[init_clu==ks]=mean(yhat_init[init_clu==ks])
  }
  init_clu_rank=cluster_reassign(y=y,new_cluster=init_clu,k=kk) # rank of initial ordinal clusters
  
  
  #########################################################################
  #########################################################################
  # iteration to find monotone clusters
  
  cluster_matrix=matrix(0,nrow(data_x),iter) # matrix to save cluster rank
  cluster_matrix_order=matrix(0,nrow(data_x),iter) # matrix to save cluster score
  fx_lst=list()
  coef_lst=list()
  y_hat_lst=list()
  stop=0 # to check the number of iteration
  # Monotone spline
  Z<-monotone.splines(data_x, num.knotsf) 
  groups <- as.vector(t(matrix(rep((1):(ncol(data_x)),(num.knotsf+2)),ncol(data_x),(num.knotsf+2))))
  objective=matrix(0,iter,1)
  
  # Initial information
  clu_rank=init_clu_rank
  y=yhat_init 
  
  for(i in 1:iter){
    #################################################################
    #################################################################
    # Fitting MS-lasso (MAM step)
    ms_lasso=MS_lasso(Xf=data_x,Zf=Z,Yf=y,num.knotsf=num.knotsf,
                           lambda=lambda)
    
    ms_coef=t(ms_lasso$coef)
    y_mean=ms_lasso$b0
    Yc=ms_lasso$Yc 
    
    # prediction
    y_hat=Z%*%ms_coef#+y_mean # estimated clusters score for all observation
    
    # Matrix of estimated monotone function of each variable
    function_mat=matrix(0,nrow(data_x),ncol(data_x))
    for(j in 1:ncol(data_x)){
      coef_j=ms_coef[groups==j,]
      z_j=Z[,groups==j]
      function_mat[,j]=z_j%*%coef_j
    }
    
    fx_lst[[i]]=function_mat
    coef_lst[[i]]=ms_lasso$coef
    y_hat_lst[[i]]=y_hat
    #################################################################
    #################################################################
    # for initial centers of UCF step
    km_mean=rep(0,kk)
    for(r in 1:kk){
      y_hat_r=y_hat[clu_rank==r,]
      km_mean[r]=mean(y_hat_r)
    }
    if(length(unique(km_mean))<kk){
      stop=stop-1
      function_mat=function_mat*0
      break
    } # If kmeans impossible because lambda is too large->stop
    
    #################################################################
    #################################################################
    # K-means
    km_new=kmeans(y_hat,centers=c(km_mean))
    o_c=km_new$cluster # cluster assignments
    clu_rank=cluster_reassign(y=y,new_cluster=o_c,k=kk) # cluster rank
    # cluster score vector
    oc_order=o_c*0 
    for(ks in 1:kk){
      oc_order[clu_rank==ks]=mean(y_hat[clu_rank==ks])
    }
    
    cluster_matrix[,i]=clu_rank # save cluster
    cluster_matrix_order[,i]=oc_order # save cluster score vector
    
    #################################################################
    #################################################################
    # Stopping rule
    if(i>1){
      if(sum(cluster_matrix[,i-1]!=cluster_matrix[,i])<(max(1,floor(delta*n))+1)){
        break
      }
    }
    #update cluster score
    y=oc_order
    # save number of iteration
    stop=stop+1
  }
  
  ########################### # retrun
  cluster_matrix=cluster_matrix[,c(1:(stop+1))]
  cluster_matrix_score=cluster_matrix_order[,c(1:(stop+1))]
  cluster=clu_rank
  
  if(i==1){
    function_mat=matrix(0,nrow(data_x),ncol(data_x))
    for(j in 1:ncol(data_x)){
      function_mat[,j]=pca_dt$rotation[j,1]*data_x[,j]
    }
    cluster=init_clu_rank
    cluster_matrix=cluster
    y_hat=PC1
  }
  
  mono=mono_check(ms_coef=ms_coef,ms_group=groups,p=ncol(data_x))
  
  if(stop>=max.iter){
    return("no convergence")
  }else if(i==1){
    return(list(iteration=0,y_hat=PC1,cluster=cluster,f_x=function_mat,Z=Z,ms_coef=ms_coef))
  }else{
    return(list(cluster=cluster, cluster_matrix=cluster_matrix,y_hat=y_hat,mono=mono,
                iteration=(stop+1),f_x=function_mat,fx_lst=fx_lst,coef_lst=coef_lst,y_hat_lst=y_hat_lst,
                y_tr_mean=y_mean,Z=Z,Yc=Yc,ms_coef=ms_coef,groups=groups))
  }
}

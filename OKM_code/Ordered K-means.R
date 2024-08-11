# Ordered K-means 
## settings :
## - First 3 variables used as ordinal variables
## - 1st and 3rd variables are related positively with cluster order
## - 2nd variable is related negatively with cluster order

## - As in paper, we use linear function as preference function


# To change the settings (ex.number of ordinal variables, direction of order),
# preference functions ("pre_l_ij", "pre_ij") should be modified.

##########################################################################
##########################  Functions  ###################################

#--------------------------------------------------------------------------
#### Preference functions

# pre_l_ij : compare preference of i-th obsevation and j-th observation for l-th variable
pre_l_ij<-function(data_i,data_j,data,l,q=0.75,direc=1){
  var_l=direc*data[,l]
  var_mm=(direc*data[,l]-min(var_l))/(max(var_l)-min(var_l))
  var_q=quantile(var_mm,0.75)
  
  pre_l_ij=direc*(c(data_i)[l]-c(data_j)[l])/(max(var_l)-min(var_l))
  
  if(pre_l_ij>var_q){
    pre_l_ij=1
  }else if(pre_l_ij<0){
    pre_l_ij=0
  }else{
    pre_l_ij=pre_l_ij/var_q
  }
  return(pre_l_ij)
}

# pre_ij : compare preference of i-th obsevation and j-th observation for based on 1~3 variables
#          using the information of settings
pre_ij<-function(data_i,data_j){
  pre_1=pre_l_ij(data_i,data_j,data=data,l=1,q=0.75,direc=1) # direc=1 -> positive
  pre_2=pre_l_ij(data_i,data_j,data=data,l=2,q=0.75,direc=-1) # direc=-1 -> negative
  pre_3=pre_l_ij(data_i,data_j,data=data,l=3,q=0.75,direc=1) # direc=1 -> positive
  
  pre_ij=mean(c(pre_1,pre_2,pre_3))
  return(pre_ij)
}

#--------------------------------------------------------------------------
####  Functions for clusteirng using preference functions

# center_order : Assign of clusters center's order based on given cluster assignment and prefernce by pre_ij
center_order<-function(data,center_mat,k){
  center_order=c(rep(0,k))
  
  for(c in 1:k){
    cent_c=center_mat[c,]
    
    for(i in 1:nrow(data)){
      d_i=data[i,]
      center_order[c]=center_order[c]+pre_ij(cent_c,d_i)-pre_ij(d_i,cent_c)
    }
  }
  return(order(center_order))
}

# okm_assign : Based on cluster center order reassinced by center_order function, assign each dataset to new clusters.
okm_assign<-function(data,center_mat,c_o){
  dist_mat=matrix(0,nrow(data),nrow(center_mat))
  
  for(i in 1:nrow(data)){
    d_i=data[i,]
    for(c in 1:nrow(center_mat)){
      cent_c=center_mat[c,]
      dist_mat[i,c_o[c]]=(pre_ij(cent_c,d_i)-pre_ij(d_i,cent_c))^2
    }
  }
  cluster=apply(dist_mat,1,which.min)
  return(cluster)
}

# center_update : updating new center based on new cluster assignment of pkm_assign
center_update<-function(data,cluster,k){
  center_mat=data[c(1:k),]*0
  
  for(c in 1:k){
    data_c=data[cluster==c,]
    dist_c_mat=matrix(0,nrow(data_c),nrow(data_c))
    
    for(i in 1:nrow(data_c)){
      for(j in 1:nrow(data_c)){
        d_i=data_c[i,]
        d_j=data_c[j,]
        
        dist_c_mat[i,j]=abs(pre_ij(d_i,d_j)-pre_ij(d_j,d_i))
      }
    }
    pre_c=apply(dist_c_mat,1,sum)
    center_mat[c,]=data_c[which.min(pre_c),]
  }
  return(center_mat)
}

#--------------------------------------------------------------------------
####  Ordered K-means 

okm_basic<-function(data,k,max.iter=500,seed=100){
  ss=seed
  data_x=as.matrix(data)
  
  # To check initial center point problem
  init_min_rate=0
  while(init_min_rate<0.05){
    set.seed(ss)
    init_sam=sample(1:nrow(data_x),k,replace=FALSE)
    center_mat=data_x[init_sam,]
    
    # initial clusters
    c_o=center_order(data_x,center_mat=center_mat,k=k)
    cluster=okm_assign(data_x,center_mat=center_mat,c_o)
    
    init_min_rate=min(table(cluster)/sum(table(cluster)))
    ss=ss+1
  }
  
  
  for(iter in 1:max.iter){
    c_o=center_order(data_x,center_mat=center_mat,k=k)
    cluster=okm_assign(data_x,center_mat=center_mat,c_o)
    center_mat_new=center_update(data_x,cluster,k)
    
    center_pre=0
    for(k_i in 1:k){
      center_pre=center_pre+pre_ij(center_mat[k_i,],center_mat_new[k_i,])
      center_pre=center_pre-pre_ij(center_mat_new[k_i,],center_mat[k_i,])
    }
    
    e_c=center_pre
    cluster_new=okm_assign(data_x,center_mat=center_mat_new,c_o)
    diff=abs(sum(cluster-cluster_new))
    
    if(e_c<0.001){
      break
    }else{
      center_mat=center_mat_new
    }
    
  }
  
  return(list(cluster=cluster,center_mat=center_mat,c_o=c_o))
  
}



##################################################################################
##########################  Clustering example  ###################################

#--------------------------------------------------------------------------
# Data Genertion
mu_or=1;n_k_lst=c(50,50,50,50);
p=20;seed=102;mu_no=5
set.seed(seed)
K=length(n_k_lst)
n=sum(n_k_lst)

for(i in 1:K){
  if(i==1){
    clu_k=c(rep(i,n_k_lst[i]))
  }else{
    clu_k=c(clu_k,rep(i,n_k_lst[i]))
  }
}

## ordinal variables
mu_1=mu_or*(2*clu_k)
mu_2=-mu_or*(exp(clu_k)/10)
mu_3=mu_or*(7*log(clu_k))

x1=mu_1+rnorm(n,0,1)
x2=mu_2+rnorm(n,0,1)
x3=mu_3+rnorm(n,0,1)

data_x=scale(cbind(x1,x2,x3))

## Data corresponding to nominal clusters
mu_no_mat=diag(K)

# Nominal cluster generation
x_no=matrix(0,n,ncol(mu_no_mat))
p_lst=runif(n,0,1)
clu_no=p_lst*0
for(i in 1:K){
  clu_no=clu_no+(p_lst>=((i-1)/K))
}
#clu_no
for(i in 1:nrow(x_no)){
  x_no[i,]=mu_no*mu_no_mat[clu_no[i],]
}

for(i in 1:ncol(x_no)){
  x_no[,i]=x_no[,i]+rnorm(nrow(x_no),0,1)
}

colnames(x_no)=paste("x_no",1:(ncol(x_no)),sep="_")

###############################
# Noise variables
p_e=p-ncol(data_x)-ncol(x_no)
noise_mat=scale(matrix(rnorm(n*p_e,0,1),n,p_e))
colnames(noise_mat)=paste("noise",1:(ncol(noise_mat)),sep="_")

data=scale(cbind(data_x,x_no,noise_mat))

#--------------------------------------------------------------------------
# example Clustering
#okm=okm_basic(data=data,k=4,seed=100)
#table(okm$cluster,clu_k) # result evaluation
#plot(okm$center_mat[,c(1:3)])


#pairs(data[,c(1,2,3)])



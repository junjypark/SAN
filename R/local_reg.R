local_reg<-function(data.mat,range=5,dis,batches,mod,standardize=TRUE){
  n=ncol(data.mat);p=nrow(data.mat)
  if(is.null(mod)){
    covariates=rep(1,n)
    q=0
  }else{
    covariates=mod[,-1]
    q=ncol(covariates)
    }
  
  batches.f=as.factor(batches)
  batch.id=unique(batches); n.batches=length(batch.id)
  count.batches=sapply(1:n.batches,function(m) sum(batches.f==batch.id[m]))
  data.std=matrix(ncol=n,nrow=p)
  data.reg=matrix(ncol=n,nrow=p)
  mod.mean=matrix(0,ncol=n,nrow=p)
  global.int<-matrix(ncol=n,nrow=p)
  sigma.mat<-matrix(ncol=n,nrow=p)
  scanner.mean<-matrix(ncol=n,nrow=p)
  res_global<-matrix(ncol=n,nrow=p)
  
  df1=sqrt((n)/(n-q-1))
  
  for (i in 1:p){
    
    y=data.mat[i,]
    global_reg<-lm(y~covariates)
    if(!is.null(mod)){
      mod.mean[i,]<-covariates%*%as.numeric(global_reg$coefficients[2:(q+1)]) 
    }
    
    global.int[i,]<- global_reg$coefficients[1]
    res_global[i,]<-residuals(global_reg)
  }
  
  for (i in 1:p){
    num_neighbors=sum(dis[i,]<=range)
    idx=order(dis[i,])[1:num_neighbors]
    
    data=res_global[idx,,drop=FALSE]*df1
    
    theta=sapply(1:n.batches, function(m) mean(data[,batches==batch.id[m],drop=FALSE]))
    scanner.mean[i,]= theta[batches.f]
    
    S=sapply(1:n.batches, function(m) sqrt((count.batches[m]-1)/count.batches[m])*sd(data[,batches==batch.id[m],drop=FALSE]))
    
    sigma.mat[i,]<-S[batches.f]
    data.reg[i,]=data[1,]-scanner.mean[i,]
    data.std[i,]=(data[1,]-scanner.mean[i,])/sigma.mat[i,]
  }
  
  if (standardize==FALSE){
    data.mat=data.reg
    sigma.mat=matrix(1,ncol=n,nrow=p)
  }else{
    data.mat=data.std
  }
  
  return(list(Xbeta=mod.mean,alpha=global.int,sigma.mat=sigma.mat,data.mat=data.mat,scanner.mean=scanner.mean))
}

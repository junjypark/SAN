local_combat<-function(data.mat,range=5,dis,batches,mod,standardize=TRUE,parallel=FALSE,ncores=2){
  n=ncol(data.mat);p=nrow(data.mat)
  if(is.null(mod)){
    q=1
  }else{
    q=ncol(mod)}
  data.std=matrix(ncol=n,nrow=p)
  mod.mean=matrix(0,ncol=n,nrow=p)
  global.int<-matrix(ncol=n,nrow=p)
  sigma.mat<-matrix(1,ncol=n,nrow=p)
  var_pooled<-matrix(ncol=n,nrow=p)
  scanner.mean<-matrix(ncol=n,nrow=p)
  scanner.mean_OLS<-matrix(ncol=n,nrow=p)
  data.combat=matrix(ncol=n,nrow=p)
  multiplicative.effect=matrix(ncol=n,nrow=p)
  multiplicative.effect_OLS=matrix(ncol=n,nrow=p)
  
  
  if (isTRUE(parallel)){
    cl = makeCluster(ncores)
    registerDoParallel(cl)
    result = foreach(i = 1:p, .packages=("neuroCombat"))%dopar%{
      num_neighbors=sum(dis[i,]<=range)
      idx=order(dis[i,])[1:num_neighbors]
      data=data.mat[idx,]
      
      combat_result<-neuroCombat(dat=data, batch=as.factor(batches), mod,verbose = FALSE)
      
      if (is.null(mod)){
        mod.mean<-rep(0,n)
      }else{
        mod.mean<-combat_result$estimates$mod.mean[1,]
      }
      
      global.int<-combat_result$estimates$stand.mean[1,]
      sigma.mat<-rep(sd(combat_result$dat.combat[1,]-mod.mean-combat_result$estimates$stand.mean[1,]),n)
      
      var_pooled<-rep(combat_result$estimates$var.pooled[1],n)
      data.std=(combat_result$dat.combat[1,]-mod.mean-combat_result$estimates$stand.mean[1,])/sigma.mat*sqrt((n-1)/(n-q))
      scanner.mean=(combat_result$estimates$gamma.star[,1])[batches]
      scanner.mean_OLS=(combat_result$estimates$gamma.hat[,1])[batches]
      multiplicative.effect=(combat_result$estimates$delta.star[,1])[batches]
      multiplicative.effect_OLS=(combat_result$estimates$delta.hat[,1])[batches]
      data.combat<-combat_result$dat.combat[1,]
      
      list(mod.mean=mod.mean,global.int=global.int,sigma.mat=sigma.mat,data.std=data.std,
           scanner.mean=scanner.mean,data.combat=data.combat, var_pooled= var_pooled,
           scanner.mean_OLS=scanner.mean_OLS,multiplicative.effect=multiplicative.effect,
           multiplicative.effect_OLS=multiplicative.effect_OLS)
    }
    
    stopCluster(cl)
    mod.mean=do.call("rbind", lapply(result, `[[`, "mod.mean"))
    global.int=do.call("rbind", lapply(result, `[[`, "global.int"))
    sigma.mat=do.call("rbind", lapply(result, `[[`, "sigma.mat"))
    data.std=do.call("rbind", lapply(result, `[[`, "data.std"))
    scanner.mean=do.call("rbind", lapply(result, `[[`, "scanner.mean"))
    data.combat=do.call("rbind", lapply(result, `[[`, "data.combat"))
    var_pooled=do.call("rbind", lapply(result, `[[`, "var_pooled"))
    scanner.mean_OLS=do.call("rbind", lapply(result, `[[`, "scanner.mean_OLS"))
    multiplicative.effect=do.call("rbind", lapply(result, `[[`, "multiplicative.effect"))
    multiplicative.effect_OLS=do.call("rbind", lapply(result, `[[`, "multiplicative.effect_OLS"))
    
    
  }else{
    
    for (i in 1:p){
      #idx=order(dis[i,])[1:neighbors]
      num_neighbors=sum(dis[i,]<=range)
      idx=order(dis[i,])[1:num_neighbors]
      data=data.mat[idx,]
      
      combat_result<-neuroCombat(dat=data, batch=as.factor(batches), mod,verbose = FALSE)
      
      if (!is.null(mod)){
        mod.mean[i,]<-combat_result$estimates$mod.mean[1,]}
      global.int[i,]<-combat_result$estimates$stand.mean[1,]
      sigma.mat[i,]<-sd(combat_result$dat.combat[1,]-mod.mean[i,]-combat_result$estimates$stand.mean[1,])
      
      var_pooled[i,]<-combat_result$estimates$var.pooled[1]
      data.std[i,]=(combat_result$dat.combat[1,]-mod.mean[i,]-combat_result$estimates$stand.mean[1,])/sigma.mat[i,]*sqrt((n-1)/(n-q))
      scanner.mean[i,]=(combat_result$estimates$gamma.star[,1])[batches]
      scanner.mean_OLS[i,]=(combat_result$estimates$gamma.hat[,1])[batches]
      multiplicative.effect[i,]=(combat_result$estimates$delta.star[,1])[batches]
      multiplicative.effect_OLS[i,]=(combat_result$estimates$delta.hat[,1])[batches]
      data.combat[i,]<-combat_result$dat.combat[1,]
    }
  }
  
  if (standardize==FALSE){
    data.mat=data.combat
    sigma.mat=matrix(1,ncol=n,nrow=p)
  }else{
    data.mat=data.std
  }
  return(list(data.combat=data.combat,Xbeta=mod.mean,alpha=global.int,sigma.mat=sigma.mat,data.mat=data.mat,scanner.mean=scanner.mean,scanner.mean_OLS=scanner.mean_OLS, multiplicative.effect= multiplicative.effect, multiplicative.effect_OLS= multiplicative.effect_OLS,var_pooled=var_pooled))
}

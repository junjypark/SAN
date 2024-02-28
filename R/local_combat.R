#' @title ComBat model for Stage 1
#' @description Helper function to perform SAN harmonization.
#' 
#' @param dat Numeric matrix with imaging features as rows,
#'                and subjects as columns.
#' @param range  Radius for Stage 1; default 5.
#' @param distMat Symmetric matrix indicating distance between imaging features.
#' @param batch Numeric or character vector specifying the batch/scanner 
#'     variable needed for harmonization. 
#' @param mod Optional model matrix for outcome of interest and 
#'     other covariates besides batch/scanner.
#' @param standardize Logical. TRUE by default; If TRUE, the output from Stage 1 is standardized.
#' @return A named \code{list} of length 10. The 1st element (\code{data.combat})
#'     contains the unstandardized ComBat harmonized data. 
#'     The 2nd element (\code{Xbeta}) contains the covariate effect. 
#'     The 3rd element (\code{alpha}) contains the global intercept.  
#'     The 4th element (\code{sigma.mat}) contains scanner-specific standard deviation. 
#'     The 5th element (\code{dat}) contains standardized ComBat harmonized data if standardize=TRUE. Otherwise, it is the same as \code{data.combat}.
#'     The 6th element (\code{scanner.mean}) contains ComBat estimated scanner-specific mean.
#'     The 7th element (\code{scanner.mean_OLS}) contains OLS estimated scanner-specific mean.
#'     The 8th element (\code{multiplicative.effect}) contains ComBat estimated multiplicative effect.
#'     The 9th element (\code{multiplicative.effect_OLS}) contains OLS estimated multiplicative effect.


local_combat=function(dat,range=5,distMat,batch,mod,standardize=TRUE){
  n=ncol(dat);p=nrow(dat)
  if(is.null(mod)){ q=1 }
  else{ q=ncol(mod) }
  data.std=matrix(ncol=n,nrow=p)
  mod.mean=matrix(0,ncol=n,nrow=p)
  global.int=matrix(ncol=n,nrow=p)
  sigma.mat=matrix(1,ncol=n,nrow=p)
  var_pooled=matrix(ncol=n,nrow=p)
  scanner.mean=matrix(ncol=n,nrow=p)
  scanner.mean_OLS=matrix(ncol=n,nrow=p)
  data.combat=matrix(ncol=n,nrow=p)
  multiplicative.effect=matrix(ncol=n,nrow=p)
  multiplicative.effect_OLS=matrix(ncol=n,nrow=p)
  
  
  for (i in 1:p){
    #idx=order(distMat[i,])[1:neighbors]
    num_neighbors=sum(distMat[i,]<=range)
    idx=order(distMat[i,])[1:num_neighbors]
    data=dat[idx,]
    
    combat_result=neuroCombat(dat=data, batch=as.factor(batch), mod,verbose = FALSE)
    
    if (!is.null(mod)){
      mod.mean[i,]=combat_result$estimates$mod.mean[1,]}
    global.int[i,]=combat_result$estimates$stand.mean[1,]
    sigma.mat[i,]=sd(combat_result$dat.combat[1,]-mod.mean[i,]-combat_result$estimates$stand.mean[1,])
    
    var_pooled[i,]=combat_result$estimates$var.pooled[1]
    data.std[i,]=(combat_result$dat.combat[1,]-mod.mean[i,]-combat_result$estimates$stand.mean[1,])/sigma.mat[i,]*sqrt((n-1)/(n-q))
    scanner.mean[i,]=(combat_result$estimates$gamma.star[,1])[batch]
    scanner.mean_OLS[i,]=(combat_result$estimates$gamma.hat[,1])[batch]
    multiplicative.effect[i,]=(combat_result$estimates$delta.star[,1])[batch]
    multiplicative.effect_OLS[i,]=(combat_result$estimates$delta.hat[,1])[batch]
    data.combat[i,]=combat_result$dat.combat[1,]
  }
  
  if (standardize==FALSE){
    dat=data.combat
    sigma.mat=matrix(1,ncol=n,nrow=p)
  }else{
    dat=data.std
  }
  return(list(data.combat=data.combat,Xbeta=mod.mean,alpha=global.int,sigma.mat=sigma.mat,dat=dat,scanner.mean=scanner.mean,scanner.mean_OLS=scanner.mean_OLS, multiplicative.effect= multiplicative.effect, multiplicative.effect_OLS= multiplicative.effect_OLS))
}

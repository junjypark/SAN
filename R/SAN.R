#' @title SAN harmonization
#' @description Main function to perform SAN harmonization.
#' 
#' @param dat Numeric matrix with imaging features as rows,
#'                and subjects as columns.
#' @param batch Numeric or character vector specifying the batch/scanner 
#'     variable needed for harmonization. 
#' @param distMat Symmetric matrix indicating distance between imaging features.
#' @param range  Radius for Stage 1; default 5.
#' @param mod Optional model matrix for outcome of interest and 
#'     other covariates besides batch/scanner.
#' @param kernel SACFs for Stage 2; default \code{kernel="mixture"}; \code{kernel="exponential"} and \code{kernel="squared exponential"} are also supported.
#' @param standardize Logical. TRUE by default; If TRUE, the output from Stage 1 is standardized.
#' @param method The integration of harmonization method; default \code{method="None"};\code{method="RELIEF"} and \code{method="CovBat"} are also supported.
#' @param parallel Logical. parallel option for ComBat in Stage 1 and conditional expectations in Stage 2; default \code{parallel=FALSE}.
#' @param ncores The number of cores when parallel computing is executed.
#' @return A named \code{list} of length 4. The 1st element (\code{data_h})
#'     contains the final harmonized data. The 2nd element (\code{data_res})
#'     contains the Stage 1 harmonized data. The 3rd element (\code{epsilon_h})
#'     contains the Stage 2 harmonized data.  The 4th element (\code{estimates}) 
#'     contains estimates and other parameters used during harmonization. 
#' 

SAN=function(dat, batch=NULL, distMat=NULL, mod=NULL, 
              range=5, kernel='mixture', standardize=TRUE, method='None',
              parallel=FALSE,ncores=2){
  
  batch.f=as.factor(batch)
  batch.id=unique(batch); n.batch=length(batch.id)
  batch.covariates=model.matrix(~batch.f-1)
  count.batch=sapply(1:n.batch,function(m) sum(batch.f==batch.id[m]))
  
  stage1_result=local_combat(dat=dat,range=range,distMat=distMat,batch=batch,mod=mod,standardize = standardize) 
  
  dat=stage1_result$dat
  
  p=nrow(dat)
  n=ncol(dat)
  
  scanner.mean=stage1_result$scanner.mean
  alpha=stage1_result$alpha
  
  sigma.mat.h=stage1_result$sigma.mat
  data.combat=stage1_result$data.combat
  
  Xbeta=stage1_result$Xbeta
  
  epsilon=list()
  for (i in 1:n.batch){ epsilon[[i]]=dat[,batch==batch.id[i]] }
  
  if (kernel=="exponential"){
    corMat.base=exp(-distMat)
    phi.hat=optimize(CovRegOptim.1,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$'minimum'
    varcomps=ObtainVarComps.1(phi.hat, epsilon=epsilon, corMat.base)
    sigma=varcomps$sigma
    tau=varcomps$tau
    harmonization_result=harmonization_step(
      kernel='exponential',sigma,tau,phi.hat,
      dat,count.batch,distMat,batch.f,batch.id,
      method=method,parallel=parallel,ncores=ncores)
    
    delta_h=harmonization_result$delta_h
    gamma_h=harmonization_result$gamma_h
    delta=harmonization_result$delta
    gamma=harmonization_result$gamma
    sigma_h=harmonization_result$sigma_h
    tau_h=harmonization_result$tau_h
    
    if (standardize){ epsilon_h=(delta_h+gamma_h)*sigma.mat.h }
    else{ epsilon_h=delta_h+gamma_h }
    
    data_h=alpha+Xbeta+epsilon_h
  }
  else if (kernel=="squared exponential"){
    corMat.base=exp(-distMat^2)
  
    phi.hat=optimize(CovRegOptim.1,interval=c(10^-5, 10),epsilon=epsilon, corMat_base=corMat.base)$'minimum'
    varcomps=ObtainVarComps.1(phi=phi.hat, epsilon=epsilon, corMat_base=corMat.base)
    sigma=varcomps$sigma
    tau=varcomps$tau
    
    harmonization_result=harmonization_step(
      kernel='squared exponential',sigma,tau,phi.hat,
      dat,count.batch,distMat,batch.f,batch.id,
      method=method,parallel=parallel,ncores=ncores)
    delta_h=harmonization_result$delta_h
    gamma_h=harmonization_result$gamma_h
    delta=harmonization_result$delta
    gamma=harmonization_result$gamma
    sigma_h=harmonization_result$sigma_h
    tau_h=harmonization_result$tau_h
    
    if (standardize){ epsilon_h=(delta_h+gamma_h)*sigma.mat.h }
    else{ epsilon_h=delta_h+gamma_h }
    
    data_h=alpha+Xbeta+epsilon_h
  }
  else if (kernel=='mixture'){
    corMat.base1=exp(-distMat)
    corMat.base2=exp(-distMat^2)
    
    phi.hat=optim(c(0.01, 0.01), CovRegOptim.2,epsilon=epsilon, corMat_base1=corMat.base1,
                  corMat_base2=corMat.base2)$par
    varcomps=ObtainVarComps.2(phi.hat, epsilon=epsilon,  corMat.base1, corMat.base2)
    sigma=matrix(varcomps$sigma,nrow=n.batch)
    tau=varcomps$tau

    harmonization_result=harmonization_step(
      kernel='mixture',sigma,tau,phi.hat,
      dat,count.batch,distMat,batch.f,batch.id,
      method=method,parallel=parallel,ncores=ncores)
    
    delta_h=harmonization_result$delta_h
    gamma_h=harmonization_result$gamma_h
    delta=harmonization_result$delta
    gamma=harmonization_result$gamma
    sigma_h=harmonization_result$sigma_h
    tau_h=harmonization_result$tau_h
    
    if (standardize){ epsilon_h=(delta_h+Reduce('+',gamma_h))*sigma.mat.h }
    else{ epsilon_h=delta_h+Reduce('+',gamma_h) }
    
    data_h=alpha+Xbeta+epsilon_h
  }
  
  estimates=list(alpha=alpha,Xbeta=Xbeta,sigma.mat.h=sigma.mat.h,phi.hat=phi.hat,sigma=sigma,tau=tau,sigma_h=sigma_h,tau_h=tau_h,delta=delta,gamma=gamma,delta_h=delta_h,gamma_h=gamma_h)
  return(list(data_h=data_h,data_res=dat,epsilon_h=epsilon_h,estimates=estimates))
}

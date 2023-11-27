harmonization_step<-function(kernel=NULL,sigma,tau,phi.hat,data.mat,count.batches,dis,batches.f,batch.id,method='None',parallel=FALSE,ncores=2){
  if (kernel=='squared exponential'){ 
    corMat.base=exp(-dis^2)}else if (kernel=="exponential"){
      corMat.base=exp(-dis)
    }else if (kernel=='mixture'){
      corMat.base1=exp(-dis)
      corMat.base2=exp(-dis^2)
    }
  
  p=nrow(data.mat)
  n=ncol(data.mat)
  n.batches=length(batch.id)
  if (kernel == "squared exponential" | kernel=="exponential"){
    delta_bat<-NULL
    gamma_bat<-NULL
    sigma_h=weighted.mean(x=sigma,w=count.batches)
    
    tau_h=weighted.mean(x=tau,w=count.batches)
    delta=delta_h=matrix(ncol=n,nrow = p)
    gamma=gamma_h=matrix(ncol=n,nrow = p)
    
    if (isTRUE(parallel)){
      cl = makeCluster(ncores)
      registerDoParallel(cl)
      results<-foreach(j = 1:n.batches) %dopar% {
        
        mat_bat=tau[j]*diag(p)+sigma[j]*corMat.base^phi.hat
        
        Sigma_inv_bat=solve(mat_bat)
        delta_bat=tau[j]*Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        gamma_bat=sigma[j]*corMat.base^phi.hat%*%Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        
        list(delta_bat=delta_bat, gamma_bat= gamma_bat)
        
      }
      
      stopCluster(cl)
      
      for (j in 1:n.batches){
        delta[,batches.f==batch.id[j]]= results[[j]]$delta_bat
        gamma[,batches.f==batch.id[j]]=  results[[j]]$gamma_bat
        
        delta_h[,batches.f==batch.id[j]]= results[[j]]$delta_bat/sqrt(tau[j])* sqrt(tau_h)
        gamma_h[,batches.f==batch.id[j]]= results[[j]]$gamma_bat/sqrt(sigma[j])* sqrt(sigma_h)
        
      } 
    }else{
      for (j in 1:n.batches){
        mat_bat=tau[j]*diag(p)+sigma[j]*corMat.base^phi.hat
        
        Sigma_inv_bat=solve(mat_bat)
        delta_bat[[j]]=tau[j]*Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        gamma_bat[[j]]=sigma[j]*corMat.base^phi.hat%*%Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        
        delta[,batches.f==batch.id[j]]= delta_bat[[j]]
        gamma[,batches.f==batch.id[j]]= gamma_bat[[j]]
        
        delta_h[,batches.f==batch.id[j]]= delta_bat[[j]]/sqrt(tau[j])* sqrt(tau_h)
        gamma_h[,batches.f==batch.id[j]]= gamma_bat[[j]]/sqrt(sigma[j])* sqrt(sigma_h)
      }
      
      
    }
    
    
    if (method=="RELIEF"){
      delta_h<-relief(dat=delta, batch=batches.f, mod =NULL,max.iter=50)$dat.relief
    }else if(method=="CovBat"){
      delta_h<-covbat(dat=delta, bat=batches.f, mod =NULL)$dat.cov
      
    }
    
  }else if (kernel=='mixture'){
    delta_bat<-NULL
    gamma1_bat<-NULL
    gamma2_bat<-NULL
    
    sigma1_h=weighted.mean(x=sigma[,1],w=count.batches)
    sigma2_h=weighted.mean(x=sigma[,2],w=count.batches)
    tau_h=weighted.mean(x=tau,w=count.batches)
    sigma_h=cbind(sigma1_h,sigma2_h)
    delta=matrix(ncol=n,nrow = p)
    gamma1=matrix(ncol=n,nrow = p)
    gamma2=matrix(ncol=n,nrow = p)
    
    delta_h=matrix(ncol=n,nrow = p)
    gamma1_h=matrix(ncol=n,nrow = p)
    gamma2_h=matrix(ncol=n,nrow = p)
    
    if (isTRUE(parallel)){
      cl = makeCluster(ncores)
      registerDoParallel(cl)
      results<-foreach(j = 1:n.batches) %dopar% {
        
        mat_bat=tau[j]*diag(p)+sigma[j,1]*exp(-phi.hat[1]*dis)+sigma[j,2]*exp(-phi.hat[2]*dis^2)
        Sigma_inv_bat=solve(mat_bat)
        delta_bat=tau[j]*Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        gamma1_bat=sigma[j,1]*exp(-phi.hat[1]*dis)%*%Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        gamma2_bat=sigma[j,2]*exp(-phi.hat[2]*dis^2)%*%Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        
        list(delta_bat=delta_bat,gamma1_bat=gamma1_bat,gamma2_bat=gamma2_bat)
      }
      stopCluster(cl)
      for (j in 1:n.batches){
        delta[,batches.f==batch.id[j]]= results[[j]]$delta_bat
        gamma1[,batches.f==batch.id[j]]= results[[j]]$gamma1_bat
        gamma2[,batches.f==batch.id[j]]= results[[j]]$gamma2_bat
        
        
        delta_h[,batches.f==batch.id[j]]= results[[j]]$delta_bat/sqrt(tau[j])* sqrt(tau_h)
        gamma1_h[,batches.f==batch.id[j]]= results[[j]]$gamma1_bat/sqrt(sigma[j,1])* sqrt(sigma1_h)
        gamma2_h[,batches.f==batch.id[j]]= results[[j]]$gamma2_bat/sqrt(sigma[j,2])* sqrt(sigma2_h)
        
      }
    }else{
      for (j in 1:n.batches){
        mat_bat=tau[j]*diag(p)+sigma[j,1]*exp(-phi.hat[1]*dis)+sigma[j,2]*exp(-phi.hat[2]*dis^2)
        Sigma_inv_bat=solve(mat_bat)
        delta_bat[[j]]=tau[j]*Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        gamma1_bat[[j]]=sigma[j,1]*exp(-phi.hat[1]*dis)%*%Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        gamma2_bat[[j]]=sigma[j,2]*exp(-phi.hat[2]*dis^2)%*%Sigma_inv_bat%*%data.mat[,batches.f==batch.id[j]]
        
        delta[,batches.f==batch.id[j]]= delta_bat[[j]]
        gamma1[,batches.f==batch.id[j]]= gamma1_bat[[j]]
        gamma2[,batches.f==batch.id[j]]= gamma2_bat[[j]]
        
        
        delta_h[,batches.f==batch.id[j]]= delta_bat[[j]]/sqrt(tau[j])* sqrt(tau_h)
        gamma1_h[,batches.f==batch.id[j]]= gamma1_bat[[j]]/sqrt(sigma[j,1])* sqrt(sigma1_h)
        gamma2_h[,batches.f==batch.id[j]]= gamma2_bat[[j]]/sqrt(sigma[j,2])* sqrt(sigma2_h)
        
      }
    }
    gamma_h<-list(gamma1_h,gamma2_h)
    gamma<-list(gamma1,gamma2)
    if (method=="RELIEF"){
      delta_h<-relief(dat=delta, batch=batches.f, mod =NULL,max.iter=50)$dat.relief
    }else if(method=="CovBat"){
      delta_h<-covbat(dat=delta, bat=batches.f, mod =NULL)$dat.cov
      
    }
    
  }
  return(list(delta=delta,delta_h=delta_h,gamma=gamma,gamma_h=gamma_h,tau_h=tau_h,sigma_h=sigma_h))
}

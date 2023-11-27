CovRegOptim.1=function(phi, epsilon, corMat_base){
  p=nrow(epsilon[[1]])
  nbatch=length(epsilon)
  nsub=sapply(epsilon, function(x) ncol(x))
  y=vector(length=2*nbatch)
  corMat=corMat_base^phi
  corMat_norm=sum(corMat^2)
  if (corMat_norm>p+1e-10){
    y1_r=sapply(epsilon, function(x) sum(x*(corMat%*%x)))
    
    comp=sapply(epsilon,function(x) sum(x^2))
    y[2*(1:nbatch)-1]=y1_r
    y[2*(1:nbatch)]=comp
    
    A=matrix(0,ncol = (2*nbatch),nrow = (2*nbatch))
    vec1=nsub*corMat_norm
    vec2=nsub*p
    for (i in 1:nbatch){
      A[(2*i-1):(2*i),(2*i-1):(2*i)]=matrix(c(vec1[i],rep(vec2[i],3)),ncol=2)
    }
    
    
    est=solve(A,y)
    sigma=est[2*(1:nbatch)-1]
    tau=est[2*(1:nbatch)]
    
    
    comp_all=cbind(comp,y1_r,sigma,nsub,tau)
    ss=sum(apply(comp_all,1,function(x) -2*(x[3]*x[2]+x[5]*x[1])))+sum(apply(comp_all, 1, function(x) x[4]*(p*x[5]*x[5]+corMat_norm*x[3]*x[3]+2*x[5]*x[3]*p)))
    
    
  } else{
    tau=-1;
    sigma=rep(-1,times=nbatch)
    ss=10^10;
  }
  return(ss)
}


ObtainVarComps.1=function(phi, epsilon, corMat_base){
  p=nrow(epsilon[[1]])
  nbatch=length(epsilon)
  nsub=sapply(epsilon, function(x) ncol(x))
  
  corMat=corMat_base^phi
  corMat_norm=sum(corMat^2)
  
  y=vector(length=2*nbatch)
  y1_r=sapply(epsilon, function(x) sum(x*(corMat%*%x)))
  
  comp=sapply(epsilon,function(x) sum(x^2))
  y[2*(1:nbatch)-1]=y1_r
  y[2*(1:nbatch)]=comp
  
  
  A=matrix(0,ncol = (2*nbatch),nrow = (2*nbatch))
  vec1=nsub*corMat_norm
  vec2=nsub*p
  for (i in 1:nbatch){
    A[(2*i-1):(2*i),(2*i-1):(2*i)]=matrix(c(vec1[i],rep(vec2[i],3)),ncol=2)
  }
  
  est=solve(A,y)
  sigma=est[2*(1:nbatch)-1]
  tau=est[2*(1:nbatch)]
  
  return(list(sigma=sigma, tau=tau))
}

###  mixture

CovRegOptim.2=function(phi, epsilon, corMat_base1, corMat_base2){
  p=nrow(epsilon[[1]])
  phi1=phi[1]
  phi2=phi[2]
  nbatch=length(epsilon)
  nsub=sapply(epsilon, function(x) ncol(x))
  y=vector(length=3*nbatch)
  corMat1=corMat_base1^phi1
  corMat2=corMat_base2^phi2
  corMat_norm1=sum(corMat1^2)
  corMat_norm2=sum(corMat2^2)
  corMat_norm12=sum(corMat1*corMat2)
  A=matrix(0,ncol = (3*nbatch),nrow = (3*nbatch))
  vec1=nsub*corMat_norm1
  vec2=nsub*corMat_norm2
  vec3=nsub*corMat_norm12
  vec4=nsub*p
  for (i in 1:nbatch){
    A[(3*i-2):(3*i),(3*i-2):(3*i)]=matrix(c(vec1[i],vec3[i],vec4[i],vec3[i],vec2[i],vec4[i],rep(vec4[i],3)),ncol=3)
  }
  if (det(A)>0){
    y1=sapply(epsilon, function(x) sum(x*(corMat1%*%x)))
    y2=sapply(epsilon, function(x) sum(x*(corMat2%*%x)))
    comp=sapply(epsilon,function(x) sum(x^2))
    y[3*(1:nbatch)-2]=y1
    y[3*(1:nbatch)-1]=y2
    y[3*(1:nbatch)]=comp
    
    
    
    est=solve(A,y)
    sigma1=est[3*(1:nbatch)-2]
    sigma2=est[3*(1:nbatch)-1]
    tau=est[3*(1:nbatch)]
    
    
    comp_all=cbind(comp,y1,y2,sigma1,sigma2,nsub,tau)
    ss=sum(apply(comp_all,1,function(x) -2*(x[4]*x[2]+x[5]*x[3]+x[7]*x[1])))+sum(apply(comp_all, 1, function(x) x[6]*(p*x[7]*x[7]+corMat_norm1*x[4]*x[4]+corMat_norm2*x[5]*x[5]+2*corMat_norm12*x[4]*x[5]+2*x[7]*x[4]*p+2*x[7]*x[5]*p)))
    
  } else{
    sigma1=rep(-1,nbatch);
    sigma2=rep(-1,nbatch);
    tau=rep(-1,nbatch);
    ss=10^10;
  }
  return(ss)
}


ObtainVarComps.2=function(phi, epsilon, corMat_base1, corMat_base2){
  p=nrow(epsilon[[1]])
  phi1=phi[1]
  phi2=phi[2]
  nbatch=length(epsilon)
  nsub=sapply(epsilon, function(x) ncol(x))
  y=vector(length=3*nbatch)
  corMat1=corMat_base1^phi1
  corMat2=corMat_base2^phi2
  corMat_norm1=sum(corMat1^2)
  corMat_norm2=sum(corMat2^2)
  corMat_norm12=sum(corMat1*corMat2)
  A=matrix(0,ncol = (3*nbatch),nrow = (3*nbatch))
  vec1=nsub*corMat_norm1
  vec2=nsub*corMat_norm2
  vec3=nsub*corMat_norm12
  vec4=nsub*p
  for (i in 1:nbatch){
    A[(3*i-2):(3*i),(3*i-2):(3*i)]=matrix(c(vec1[i],vec3[i],vec4[i],vec3[i],vec2[i],vec4[i],rep(vec4[i],3)),ncol=3)
  }
  
  y1=sapply(epsilon, function(x) sum(x*(corMat1%*%x)))
  y2=sapply(epsilon, function(x) sum(x*(corMat2%*%x)))
  comp=sapply(epsilon,function(x) sum(x^2))
  y[3*(1:nbatch)-2]=y1
  y[3*(1:nbatch)-1]=y2
  y[3*(1:nbatch)]=comp
  
  est=solve(A,y)
  sigma1=est[3*(1:nbatch)-2]
  sigma2=est[3*(1:nbatch)-1]
  tau=est[3*(1:nbatch)]
  
  
  return(list(sigma=c(sigma1,sigma2), tau=tau))
}

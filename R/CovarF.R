#' @title CovarF statistic
#' @description a function for computing CovarF statistics of the brain map
#' @param data.mat Numeric matrix with imaging features as rows,
#'  and subjects as columns.
#' @param range  Radius for the local neighbors; default 5.
#' @param dis Symmetric matrix indicating distance between imaging features.
#' @param batch Numeric or character vector specifying the batch/scanner 
#'     variable needed for harmonization. 
#' @param mod Optional model matrix for outcome of interest and 
#'     other covariates besides batch/scanner.
#' @param cov_true Positive definite matrix for shared covariances across batches/scanners.
#' @return A vector of CovarF statistics.

CovarF<-function(data.mat, batch, distMat, mod, range=5,cov_true){
  n=ncol(data.mat)
  p=nrow(data.mat)
  F_stat=vector(length = p)
  batch.f=as.factor(batch)
  batch.id=unique(batch)
  n.batch=length(batch.id)
  count.batch=sapply(1:n.batch,function(m) sum(batch.f==batch.id[m]))
  
  for (i in 1:p){
    num_neighbors=sum(dis[i,]<=range)
    idx=order(dis[i,])[1:num_neighbors]
    data=data.mat[idx,]
    cov_local=cov_true[idx,idx]
    res_data<-apply(t(data),2,function(x) lm(x~mod+batch.f)$residuals)
    
    res_data.ls=lapply(1:n.batch,function(i) res_data[batch==batch.id[i],])
    cov_data.ls=lapply(1:n.batch,function(i) cov(res_data.ls[[i]]))
    
    data_covdiff=sapply(1:n.batch,function(i) sum(count.batch[i]*(cov_data.ls[[i]]-cov_local)^2))
    
    cov_sub<-NULL
    for (j in 1:n.batch){
      cov_sub[[j]]=sapply(1:count.batch[j],function(i) outer(res_data.ls[[j]][i,],res_data.ls[[j]][i,],FUN = '*'))
    }
    
    SSE=sum(sapply(1:n.batch,function(i) sum((cov_sub[[i]]-c(cov_data.ls[[i]]))^2)))
    SSB=sum(data_covdiff)
    F_stat[i]=SSB/(SSE/(sum(count.batch)-2))
  }
  return(F_stat)
}


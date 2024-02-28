\name{CASH score}
\alias{CASH}
\title{CASH score}
\usage{
CASH(
data.mat,
batch, 
distMat,
mod, 
range=5,
cov_true=NULL)
}

\arguments{
\item{data.mat}{Numeric matrix with imaging features as rows, and subjects as columns.}

\item{range}{Radius for the local neighbors; default 5.}

\item{distMat}{Symmetric matrix indicating distance between imaging features.}

\item{batch}{Numeric or character vector specifying the batch/scanner variable needed for harmonization.}

\item{mod}{Optional model matrix for outcome of interest and other covariates besides batch/scanner.}

\item{cov_true}{Optional shared covariance matrix across batches/scanners. If not specified, this covariance matrix is estimated from pooled \code{data.mat}.
}



}

\value{  A vector of CASH scores.
}

\description{
A function for computing CASH scores of the brain map.
}

\examples{
CASH_score = CASH(data.mat=data.mat, batch=batch, distMat=distMat, mod=mod)

}

\references{
Zhang, R., Chen, L., Oliver, L. D., Voineskos, A. N., Park, J. Y. (2023+). SAN: mitigating spatial covariance heterogeneity in cortical thickness data collected from multiple scanners or sites. }


\name{SAN}
\alias{SAN}
\title{SAN harmonization}
\usage{
SAN(
  dat, 
  batch=NULL, 
  distMat=NULL, 
  mod=NULL, 
  range=5, 
  kernel='mixture', 
  standardize=TRUE, 
  method='None',
  parallel=FALSE,
  ncores=2
)
}

\arguments{
\item{dat}{Numeric matrix with imaging features as rows, and subjects as columns.}

\item{batch}{Numeric or character vector specifying the batch/scanner variable needed for harmonization.}

\item{distMat}{Symmetric matrix indicating distance between imaging features.}

\item{range}{Radius for Stage 1; default 5.}

\item{mod}{Optional model matrix for outcome of interest and other covariates besides batch/scanner.}

\item{kernel}{SACFs for Stage 2; default \code{kernel="mixture"}; \code{kernel="exponential"} and \code{kernel="squared exponential"} are also supported.}

\item{standardize}{Logical. TRUE by default; If TRUE, the output from Stage 1 is standardized.}

\item{method}{The integration of harmonization method; default \code{method="None"};\code{method="RELIEF"} and \code{method="CovBat"} are also supported.}


\item{parallel}{Logical. parallel option for ComBat in Stage 1 and conditional expectations in Stage 2; default \code{parallel=FALSE}.}

\item{ncores}{The number of cores when parallel computing is executed.}

}

\value{ A list of the following is returned. 
\item{data_h}{contains the final harmonized data.} \item{data_res}{contains the Stage 1 harmonized data.} 
\item{epsilon_h}{contains the Stage 2 harmonized data..} \item{estimates}{contains estimates and other parameters used during harmonization.} 
}

\description{
Main function to perform SAN harmonization.
}

\examples{
fit = SAN(dat=dat, batch=batch, mod=mod, distMat=distMat)

}

\references{
Zhang, R., Chen, L., Oliver, L. D., Voineskos, A. N., Park, J. Y. (2023+). SAN: mitigating spatial covariance heterogeneity in cortical thickness data collected from multiple scanners or sites. }


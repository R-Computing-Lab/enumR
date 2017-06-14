#' @title Generate Effects matrix
#'
#' @description
#' \code{GenStructure} returns an effects matrix, based on intercorrelations and number of factors
#'
#'
#' @param nfactors Number of factors
#' @param rfactors Intercorrelations between factors. Default is 0
#' @param r_norm If TRUE, generate normally distributed factor loadings with  mean \code{loading} and sd \code{r_norm_sd}. Does not check if correlation matrix is positive definite.
#' @param r_norm_sd If \code{r_norm} is TRUE, standard deviation for normally distributed intercorrelations.
#' @return \code{matrix} of factors and their intercorrelations.


GenStructure = function(nfactors=5,rfactors=0,r_norm=FALSE,r_norm_sd=.015){
  if(!r_norm){
  Structure =matrix(rep(rfactors,nfactors^2),ncol=nfactors)
  }else{Structure =diag(nfactors)
                         Structure[upper.tri(Structure, diag=FALSE)] <- rnorm(mean=rfactors,sd=r_norm_sd,n=(nfactors^2-nfactors)/2)
                         Structure <- round(Structure + t(Structure) - diag(diag(Structure)) ,digits=2)}
  return(Structure)
}

#' @title Generate Factor Loadings Matrix
#'
#' @description
#' \code{GenFactorMatrix} returns a matrix of factor loadings, based on number of items, factors, and loadings
#'
#'
#' @param items vector of total number of items per factor
#' @param loading Factor loading magnitude. Default is .5
#' @param loading_norm If TRUE, generate factor loadings that average to \code{loading}
#' @param loading_norm_sd If \code{loading_norm} is TRUE, standard deviation of loadings.
#' @param nfactors Number of factors
#' @param itemsR Number of items per factor that are reverse scored
#' @return \code{matrix} of n factors by i items

GenFactorMatrix = function(nfactors = 5, items= c(5,5,5,5,5), itemsR=c(2,2,2,2,2), loading=.5,loading_norm=FALSE,loading_norm_sd=.025){ # Generate Factor Loadings

  Loadings = matrix(rep(0,sum(items)*nfactors),
  ncol=nfactors,byrow=TRUE)
  if(!loading_norm){
    for(i in 1:nfactors){
    Loadings[(sum(items[1:i-1])+1):sum(items[1:i]),i]<-rep(loading,items[i])

    if(itemsR[i]>0){ #reverse score items
    Loadings[(sum(items[1:i-1])+1):(sum(items[1:i-1])+itemsR[i]),i]<- -rep(loading,itemsR[i])}
    }
    } else{
      for(i in 1:nfactors){
        Loadings[(sum(items[1:i-1])+1):sum(items[1:i]),i]<-round(rnorm(items[i], mean = loading,sd=loading_norm_sd),digits=2)

        if(itemsR[i]>0){ #reverse score items
          Loadings[(sum(items[1:i-1])+1):(sum(items[1:i-1])+itemsR[i]),i]<- - Loadings[(sum(items[1:i-1])+1):(sum(items[1:i-1])+itemsR[i]),i]
    }
      }
    }
  return(Loadings)
}


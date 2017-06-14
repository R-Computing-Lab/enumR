#' @title Generate Data from factor loadings
#'
#' @description
#' \code{GenData} returns a data.frame of observed scores based on matrix of factor loadings, effects matrix. The heart of this function was adapted from \url{http://personality-project.org/r/r.datageneration.html}

#'
#' @details
#'
#'
#' @param fmodel Matrix of factor loadings for all items
#' @param effect effects matrix
#' @param samplesize number of observations generated
#' @param names names of variables
#' @return data.frame of observed scores
#'

GenData = function(fmodel,effect,samplesize,inames=NULL) {   # define a general function in terms of a factor model and an effects matrix

  numberofvariables = dim(fmodel)[1]        #problem size determined by input to the function
  if(is.null(inames)){inames=paste0("item_",1:numberofvariables)}
  numberoflatent = dim(fmodel)[2]
  tmodel = t(fmodel)      #transpose of model
  communality=diag(fmodel%*%tmodel)       #find how much to weight true scores and errors given the measurement model
  uniqueness=1-communality
  errorweight=diag(sqrt(uniqueness))      #how much to weight the errors

  latentscores = matrix(rnorm(samplesize*(numberoflatent)),samplesize) #create true scores for the latent variables
  latentscores = latentscores%*%effect
  truescores = latentscores%*%tmodel
  error = matrix(rnorm(samplesize*(numberofvariables)),samplesize)  #create normal error scores
  error = error%*%errorweight
  observedscore = truescores+error
  observedscore = data.frame(observedscore)
  names(observedscore)<-inames
  return(observedscore) }       #end of function mes


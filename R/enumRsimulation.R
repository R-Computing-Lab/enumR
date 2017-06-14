#' @title enumR simulation function
#'
#' @description This function simulates factor analyses based on numerous parameters, listed below.
#'
#' @param cor	What kind of correlation to find, defaults to covariance matrix, but see fa for the choices
#' @param effectmatrix  OPTIONAL matrix of item loadings by factor. If not provided, one is generated from \code{items} and \code{nfactors}.
#' @param items_p_f  number of items per factor
#' @param itemsR_p_f number of reverse-scored items per factor
#' @param custom_item If TRUE, overrides \code{items_p_f} and \code{itemsR_p_f}. Primary used for importing from shiny interface.
#'  @param items OPTIONAL vector of total number of items per factor
#' @param loading Factor loading magnitude. Default is .5
#' @param loading_norm If TRUE, generate factor loadings that average to \code{loading}
#' @param loading_norm_sd If \code{loading_norm} is TRUE, standard deviation of loadings.
#' @param method factoring method – fm="pa" Principal Axis Factor Analysis, fm = "minres" minimum residual (OLS) factoring fm="ml" Maximum Likelihood FA, fm="pc" Principal Components"
#' @param samplesize Number of observations
#' @param ndatasets Number of datasets generated
#' @param nfactors Number of factors to extract.
#' @param patternmatrix OPTIONAL matrix of factor loadings. If not provided, one is generated from \code{rfactors} and \code{nfactors}.
#' @param rfactors Intercorrelations between factors. Default is 0
#' @param r_norm If TRUE, generate normally distributed factor loadings with  mean \code{loading} and sd \code{r_norm_sd}. Does not check if correlation matrix is positive definite.
#' @param r_norm_sd If \code{r_norm} is TRUE, standard deviation for normally distributed intercorrelations.
#' @param itemsR Number of items per factor that are reverse scored
#' @param rotation what rotation to use c("none", "varimax", "oblimin","promax")
#' @param seed numeric value for setting a seed. Allows results to be reproduced.
#' @param true_model If TRUE, Estimates model using \code{nfactors} only. Is used in \code{enumRsimulation} function. Default is FALSE.
#' @param use	If doing covariances or Pearson R, should we use "pairwise" or "complete cases"
#' @param  f1_items Total number of items on factor 1
#' @param  f2_items Total number of items on factor 2
#' @param  f3_items Total number of items on factor 3
#' @param  f4_items Total number of items on factor 4
#' @param  f5_items Total number of items on factor 5
#' @param  f6_items Total number of items on factor 6
#' @param  f7_items Total number of items on factor 7
#' @param  f8_items Total number of items on factor 8
#' @param  f9_items Total number of items on factor 9
#' @param  f10_items Total number of items on factor 10
#' @param  f1_itemsR Number of reverse-scored items on factor 1
#' @param  f2_itemsR Number of reverse-scored items on factor 2
#' @param  f3_itemsR Number of reverse-scored items on factor 3
#' @param  f4_itemsR Number of reverse-scored items on factor 4
#' @param  f5_itemsR Number of reverse-scored items on factor 5
#' @param  f6_itemsR Number of reverse-scored items on factor 6
#' @param  f7_itemsR Number of reverse-scored items on factor 7
#' @param  f8_itemsR Number of reverse-scored items on factor 8
#' @param  f9_itemsR Number of reverse-scored items on factor 9
#' @param  f10_itemsR Number of reverse-scored items on factor 10
#' ...	parameters to pass to the factor analysis program The most important of these is if using a correlation matrix is covmat= xx
#' @return dataframe

enumRsimulation<-function(seed=12345,
                ndatasets=200,
                patternmatrix=NULL,
                effectmatrix=NULL,
                nfactors =5,
                loading=.5,
                items=NULL,
                items_p_f=5,
                itemsR_p_f=2,
                itemsR=NULL,
                loading_norm=FALSE,
                loading_norm_sd=.05,
                # effects matrix
                rfactors=0,
                r_norm =FALSE,
                r_norm_SD =.015,
                samplesize=300,
                method="ml",
                rotation="oblimin",
                custom_item=FALSE,
                f1_items=NULL,
                f2_items = NULL,
                f3_items = NULL,
                f4_items = NULL,
                f5_items = NULL,
                f6_items = NULL,
                f7_items = NULL,
                f8_items = NULL,
                f9_items = NULL,
                f10_items = NULL,
                f1_itemsR= NULL,
                f2_itemsR = NULL,
                f3_itemsR = NULL,
                f4_itemsR = NULL,
                f5_itemsR = NULL,
                f6_itemsR = NULL,
                f7_itemsR = NULL,
                f8_itemsR = NULL,
                f9_itemsR = NULL,
                f10_itemsR = NULL
                ){
###prep
set.seed(seed)
require(psych)
require(GPArotation)
require(GAIPE)

simulation<-ndatasets

# cleanup items
if(is.null(items)){
if(custom_item){
  items=c(f1_items,
          f2_items,
          f3_items,
          f4_items,
          f5_items,
          f6_items,
          f7_items,
          f8_items,
          f9_items,
          f10_items)

  itemsR=c(f1_itemsR,
          f2_itemsR,
          f3_itemsR,
          f4_itemsR,
          f5_itemsR,
          f6_itemsR,
          f7_itemsR,
          f8_itemsR,
          f9_itemsR,
          f10_itemsR)

}else{
  items=rep(items_p_f,nfactors)
  itemsR=rep(itemsR_p_f,nfactors)
}}





## Results
meganames<-c( "dof","chisq","prob","sqresid","fit", "RMSEA","RMSEA_lower","RMSEA_upper","BIC","SABIC","null.model","null.dof","null.chisq","complex","eChisq","SRMR","eCRMS","eBIC","TLI","fa.fit","fa.fit.off","objective","ENull","eSABIC","null.chisq1" , "TLI.m","NFI.m","CFI.m","RMSEA_lower.m","RMSEA_upper.m","eRMS","cfit.1","cfit.2","cfit.3","cfit.4","cfit.5","cfit.6","cfit.7","cfit.8","cfit.9","cresidual.1","cresidual.2","cresidual.3","cresidual.4","cresidual.5","cresidual.6","cresidual.7","cresidual.8","cresidual.9","map.values","rfactors","samplesize","loadings","dataset","error","factor")
mega<-data.frame(matrix(rep(0,length(meganames)),nrow=1))
names(mega)=meganames
mega_null<-mega
mega_null$error<- 1

percentiles<-data.frame(c(1))
names(percentiles)<-"ds"


if(is.null(patternmatrix)){
  fmodel=GenFactorMatrix(nfactors=nfactors, items= items[1:nfactors], itemsR=itemsR[1:nfactors], loading=loading,loading_norm=loading_norm,loading_norm_sd=loading_norm_sd)
}else{fmodel=patternmatrix}
  if(is.null(effectmatrix)){
    smodel=GenStructure(nfactors=nfactors,rfactors=rfactors,r_norm=r_norm,r_norm_sd=r_norm_sd)
    }else{smodel=effectmatrix}
## Gen Data
for(d in 1:simulation) {
        data=GenData(fmodel=fmodel,effect=smodel,samplesize=samplesize)
        cov<-cov(data)
        test<-	try(results<- enumR(cov,n=nfactors,samplesize=samplesize,covar=TRUE,truemodel = TRUE,fm = method), silent=TRUE)
        if ('try-error' %in% class(test)) {mega[d,]<- mega_null
        } else {mega[d,]<- results
        mega$factor[d]<-nfactors
        mega$error[d]<- 0
        }
mega$dataset[d]<-d

        print(d)
}
mega$rfactors<-rfactors
mega$samplesize<-samplesize
mega$loadings<-loading

ds_true<-mega[mega$factor==nfactors,]

percentiles$RMSEA_.95<-quantile(ds_true$RMSEA,.95)
percentiles$MAP_.95<-quantile(ds_true$map.values,.95)
#percentiles$chisq_.95<-quantile(ds_true$chisq,.95)
#percentiles$BIC_.95<-quantile(ds_true$BIC,.95)
percentiles$SRMR_.95<-quantile(ds_true$SRMR,.95)
#percentiles$SABIC_.95<-quantile(ds_true$SABIC,.95)
#percentiles$TLI_.95<-quantile(ds_true$TLI,.95)
percentiles$chisq_sig_.95<-quantile(ds_true$prob,.95)
#percentiles$RMSEA_.05<-quantile(ds_true$RMSEA,.05)
#percentiles$MAP_.05<-quantile(ds_true$map.values,.05)
#percentiles$chisq_.05<-quantile(ds_true$chisq,.05)
#percentiles$BIC_.05<-quantile(ds_true$BIC,.05)
#percentiles$SRMR_.05<-quantile(ds_true$SRMR,.05)
#percentiles$SABIC_.05<-quantile(ds_true$SABIC,.05)
percentiles$TLI_.05<-quantile(ds_true$TLI,.05)
percentiles$chisq_sig_.05<-quantile(ds_true$prob,.05)

return(list(percentiles,mega))
}

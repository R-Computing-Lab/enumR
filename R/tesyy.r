library(roxygen2)
ffd=GenFactorMatrix()
effectmatrix=GenStructure()

df=GenData(ffd,effectmatrix,samplesize=1000)

enumR(cov(df),nfactors=5,truemodel = TRUE )

(maybe=enumRsimulation(samplesize=500, loading_norm=FALSE,
                       loading_norm_sd=.05,
                       r_norm =FALSE,
                       r_norm_SD =.015))[1]

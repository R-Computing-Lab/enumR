#' @title enumR function
#'
#' @description The enumR function is a heavily modified version of the \code{vss} function from \code{Psych} Package]. Both functions are designed to help researchers determine the number of factors in their exploratory factor analysis. This function includes a variarity of methods.
#'
#'
#' @param cor	What kind of correlation to find, defaults to covariance matrix, but see fa for the choices
#' @param diagonal	Should we fit the diagonal as well
#' @param method	 factoring method – fm="pa" Principal Axis Factor Analysis, fm = "minres" minimum residual (OLS) factoring fm="ml" Maximum Likelihood FA, fm="pc" Principal Components"
#' @param samplesize	 Number of observations if doing a factor analysis of correlation matrix. This value is ignored by VSS but is necessary for the ML factor analysis packag
#' @param nfactors Maximum number of factors to extract. Needs to be larger than expected.
#' @param rotation what rotation to use c("none", "varimax", "oblimin","promax")
#' @param true_model If TRUE, Estimates model using \code{nfactors} only. Is used in \code{enumRsimulation} function. Default is FALSE.
#' @param use	If doing covariances or Pearson R, should we use "pairwise" or "complete cases"
#' @param x a correlation matrix or a data matrix
#' ...	parameters to pass to the factor analysis program The most important of these is if using a correlation matrix is covmat= xx
#' @return dataframe


enumR<-function(x, nfactors = 8, rotation = "oblimin", diagonal = FALSE, method = "ml",samplesize=NULL,use="pairwise",cor="cov",truemodel=FALSE,...){
  require(psych)
  require(GAIPE)
  n=nfactors; fm=method
 # cl <- match.call()
  if (rotation == "oblimin") {
    if (!requireNamespace("GPArotation")) {
      stop("You must have GPArotation installed to use oblimin rotation")
    }
  }

  old_rotation = rotation

  if(truemodel){
    complexrow <- function(x, c) {
      n = length(x)
      temp <- x
      x <- rep(0, n)
      for (j in 1:c) {
        locmax <- which.max(abs(temp))
        x[locmax] <- sign(temp[locmax]) * max(abs(temp))
        temp[locmax] <- 0
      }
      return(x)
    }
    complexmat <- function(x, c) {
      nrows <- dim(x)[1]
      ncols <- dim(x)[2]
      for (i in 1:nrows) {
        x[i, ] <- complexrow(x[i, ], c)
      }
      return(x)
    }
    map <- function(x, n) {
      nvar <- dim(x)[2]
      min.partial <- rep(NA, n)
      e <- eigen(x)
      evect <- e$vectors
      comp <- evect %*% diag(sqrt(e$values))
      if (n >= nvar) {
        n1 <- nvar - 1
      } else {
        n1 <- n
      }

        c11.star <- x - comp[, 1:n1] %*% t(comp[, 1:n1])
        d <- diag(1/sqrt(diag(c11.star)))
        rstar <- d %*% c11.star %*% d
        diag(rstar) <- 0
        min.partial <- sum(rstar * rstar)/(nvar * (nvar -1))

      return(min.partial)
    }
    if (dim(x)[2] < n)
      n <- dim(x)[2]
    complexfit <- array(0, dim = c(1, 1))
    complexresid <- array(0, dim = c(1, 1))
    vss.df <- data.frame(dof = 0, chisq = NA, prob = NA,sqresid = NA, fit = NA, RMSEA = NA,RMSEA_lower=NA,RMSEA_upper=NA, BIC = NA, SABIC = NA, null.model = NA, null.dof = NA, null.chisq = NA,complex = NA, eChisq = NA, SRMR = NA, eCRMS = NA, eBIC = NA, TLI=NA,fa.fit=NA,fa.fit.off=NA, objective=NA,ENull=NA,eSABIC=NA,null.chisq1=NA,TLI.m=NA, NFI.m=NA, CFI.m=NA,RMSEA_lower.m=NA,RMSEA_upper.m=NA)

    if (!is.matrix(x))
      x <- as.matrix(x)

    if (is.null(samplesize)) {
      message("samplesize was not specified and was arbitrarily set to 1000.  This only affects the chi square values.")
      samplesize <- 1000}

    map.values <- map(x, n)
    if (n > dim(x)[2]) {
      n <- dim(x)[2]}


      PHI <- diag(n)
      if (n < 2) {
        (rotation = "none")
      }
      else {
        rotation = old_rotation
      }

      f <- fa(x, n, rotate = rotation, n.obs = samplesize, warnings = FALSE,
              fm = fm, scores = "none", cor = cor, ...)

        original <- x
        sqoriginal <- original * original
        totaloriginal <- sum(sqoriginal) - diagonal *
          sum(diag(sqoriginal))


      load <- as.matrix(f$loadings)
      model <- load %*% PHI %*% t(load)
      residual <- original - model
      sqresid <- residual * residual
      totalresid <- sum(sqresid) - diagonal * sum(diag(sqresid))
      fit <- 1 - totalresid/totaloriginal

      vss.df[ "dof"] <- f$dof
      vss.df[ "chisq"] <- f$STATISTIC
      vss.df[ "prob"] <- f$PVAL
      vss.df["eChisq"] <- f$chi
      vss.df["ENull"] <- f$ENull
      vss.df["SRMR"] <- f$rms
      vss.df["eRMS"] <- f$rms
      vss.df["eCRMS"] <- f$crms
      vss.df["eBIC"] <- f$EBIC
      vss.df["eSABIC"] <- f$ESABIC

      vss.df["TLI"]<-f$TLI
      vss.df["TLI.m"]<-(f$null.chisq/f$null.dof - f$STATISTIC/f$dof)/(f$null.chisq/f$null.dof - 1)

      vss.df["NFI.m"]<-(f$null.chisq - f$STATISTIC)/(f$null.chisq)

      vss.df["CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)

      vss.df["CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)

      vss.df["null.model"] = f$null.model
      vss.df["null.dof"] = f$null.dof

      vss.df["fa.fit"] = f$fit
      vss.df["fa.fit.off"] = f$fit.off
      vss.df["objective"] = f$objective
      if (!is.null(f$chisq)) {
        vss.df["null.chisq1"] = f$chisq
        vss.df["null.chisq"] = f$null.chisq
      }
      if (!is.null(f$RMSEA[1])) {
        vss.df[ "RMSEA"] <- f$RMSEA[1]
        vss.df[ "RMSEA_upper"] <- f$RMSEA[3]
        vss.df[ "RMSEA_lower"] <- f$RMSEA[2]

        RMSEA.m<-CI.RMSEA(f$RMSEA[1],N=samplesize,df=f$dof,clevel=1-f$RMSEA[4])

        vss.df[ "RMSEA_upper.m"] <- RMSEA.m$Upper.CI
        vss.df[ "RMSEA_lower.m"] <- RMSEA.m$Lower.CI

      }
      else {
        vss.df[ "RMSEA"] <- NA
      }

      if (!is.null(f$BIC)) {
        vss.df[ "BIC"] <- f$BIC
      }
      else {
        vss.df[ "BIC"] <- NA
      }
      if (!is.null(f$SABIC)) {
        vss.df[ "SABIC"] <- f$SABIC
      }
      else {
        vss.df[ "SABIC"] <- NA
      }
      if (!is.null(f$complexity)) {
        vss.df[ "complex"] <- mean(f$complexity)
      }
      else {
        vss.df[ "complex"] <- NA
      }

      vss.df[ "sqresid"] <- totalresid
      vss.df[ "fit"] <- fit
      for (c in 1:n) {
        simpleload <- complexmat(load, c)
        model <- simpleload %*% PHI %*% t(simpleload)
        residual <- original - model
        sqresid <- residual * residual
        totalsimple <- sum(sqresid) - diagonal * sum(diag(sqresid))
        simplefit <- 1 - totalsimple/totaloriginal
        complexresid[c] <- totalsimple
        complexfit[c] <- simplefit
      }

  }else{
  complexrow <- function(x, c) {
    n = length(x)
    temp <- x
    x <- rep(0, n)
    for (j in 1:c) {
      locmax <- which.max(abs(temp))
      x[locmax] <- sign(temp[locmax]) * max(abs(temp))
      temp[locmax] <- 0
    }
    return(x)
  }
  complexmat <- function(x, c) {
    nrows <- dim(x)[1]
    ncols <- dim(x)[2]
    for (i in 1:nrows) {
      x[i, ] <- complexrow(x[i, ], c)
    }
    return(x)
  }
  map <- function(x, n) {
    nvar <- dim(x)[2]
    min.partial <- rep(NA, n)
    e <- eigen(x)
    evect <- e$vectors
    comp <- evect %*% diag(sqrt(e$values))
    if (n >= nvar) {
      n1 <- nvar - 1
    } else {
      n1 <- n
    }
    for (i in 1:n1) {
      c11.star <- x - comp[, 1:i] %*% t(comp[, 1:i])
      d <- diag(1/sqrt(diag(c11.star)))
      rstar <- d %*% c11.star %*% d
      diag(rstar) <- 0
      min.partial[i] <- sum(rstar * rstar)/(nvar * (nvar -1))
    }
    return(min.partial)
  }
  if (dim(x)[2] < n)
    n <- dim(x)[2]
  complexfit <- array(0, dim = c(n, n))
  complexresid <- array(0, dim = c(n, n))
  vss.df <- data.frame(dof = rep(0, n), chisq = NA, prob = NA,sqresid = NA, fit = NA, RMSEA = NA,RMSEA_lower=NA,RMSEA_upper=NA, BIC = NA, SABIC = NA, null.model = NA, null.dof = NA, null.chisq = NA,complex = NA, eChisq = NA, SRMR = NA, eCRMS = NA, eBIC = NA, TLI=NA,fa.fit=NA,fa.fit.off=NA, objective=NA,ENull=NA,eSABIC=NA,null.chisq1=NA,TLI.m=NA, NFI.m=NA, CFI.m=NA,RMSEA_lower.m=NA,RMSEA_upper.m=NA)

  if (!is.matrix(x))
    x <- as.matrix(x)

  if (is.null(samplesize)) {
    message("samplesize was not specified and was arbitrarily set to 1000.  This only affects the chi square values.")
    samplesize <- 1000}

  map.values <- map(x, n)
  if (n > dim(x)[2]) {
    n <- dim(x)[2]}

  for (i in 1:n) {
    PHI <- diag(i)
    if (i < 2) {
      (rotation = "none")
    }
    else {
      rotation = old_rotation
    }

    f <- fa(x, i, rotate = rotation, n.obs = samplesize, warnings = FALSE,
            fm = fm, scores = "none", cor = cor, ...)
    if (i == 1) {
      original <- x
      sqoriginal <- original * original
      totaloriginal <- sum(sqoriginal) - diagonal *
        sum(diag(sqoriginal))
    }

    load <- as.matrix(f$loadings)
    model <- load %*% PHI %*% t(load)
    residual <- original - model
    sqresid <- residual * residual
    totalresid <- sum(sqresid) - diagonal * sum(diag(sqresid))
    fit <- 1 - totalresid/totaloriginal

    vss.df[i, "dof"] <- f$dof
    vss.df[i, "chisq"] <- f$STATISTIC
    vss.df[i, "prob"] <- f$PVAL
    vss.df[i, "eChisq"] <- f$chi
    vss.df[i, "ENull"] <- f$ENull
    vss.df[i, "SRMR"] <- f$rms
    vss.df[i, "eRMS"] <- f$rms
    vss.df[i, "eCRMS"] <- f$crms
    vss.df[i, "eBIC"] <- f$EBIC
    vss.df[i, "eSABIC"] <- f$ESABIC

    vss.df[i,"TLI"]<-f$TLI
    vss.df[i,"TLI.m"]<-(f$null.chisq/f$null.dof - f$STATISTIC/f$dof)/(f$null.chisq/f$null.dof - 1)

    vss.df[i,"NFI.m"]<-(f$null.chisq - f$STATISTIC)/(f$null.chisq)

    vss.df[i,"CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)

    vss.df[i,"CFI.m"]<-(f$null.chisq-f$null.dof - f$STATISTIC-f$dof)/(f$null.chisq-f$null.dof)

    vss.df[i,"null.model"] = f$null.model
    vss.df[i,"null.dof"] = f$null.dof

    vss.df[i,"fa.fit"] = f$fit
    vss.df[i,"fa.fit.off"] = f$fit.off
    vss.df[i,"objective"] = f$objective
    if (!is.null(f$chisq)) {
      vss.df[i,"null.chisq1"] = f$chisq
      vss.df[i,"null.chisq"] = f$null.chisq
    }
    if (!is.null(f$RMSEA[1])) {
      vss.df[i, "RMSEA"] <- f$RMSEA[1]
      vss.df[i, "RMSEA_upper"] <- f$RMSEA[3]
      vss.df[i, "RMSEA_lower"] <- f$RMSEA[2]

      RMSEA.m<-CI.RMSEA(f$RMSEA[1],N=samplesize,df=f$dof,clevel=1-f$RMSEA[4])

      vss.df[i, "RMSEA_upper.m"] <- RMSEA.m$Upper.CI
      vss.df[i, "RMSEA_lower.m"] <- RMSEA.m$Lower.CI

    }
    else {
      vss.df[i, "RMSEA"] <- NA
    }

    if (!is.null(f$BIC)) {
      vss.df[i, "BIC"] <- f$BIC
    }
    else {
      vss.df[i, "BIC"] <- NA
    }
    if (!is.null(f$SABIC)) {
      vss.df[i, "SABIC"] <- f$SABIC
    }
    else {
      vss.df[i, "SABIC"] <- NA
    }
    if (!is.null(f$complexity)) {
      vss.df[i, "complex"] <- mean(f$complexity)
    }
    else {
      vss.df[i, "complex"] <- NA
    }

    vss.df[i, "sqresid"] <- totalresid
    vss.df[i, "fit"] <- fit
    for (c in 1:i) {
      simpleload <- complexmat(load, c)
      model <- simpleload %*% PHI %*% t(simpleload)
      residual <- original - model
      sqresid <- residual * residual
      totalsimple <- sum(sqresid) - diagonal * sum(diag(sqresid))
      simplefit <- 1 - totalsimple/totaloriginal
      complexresid[i, c] <- totalsimple
      complexfit[i, c] <- simplefit
    }
  }
  }
  vss.stats <- data.frame(vss.df, cfit = complexfit, cresidual = complexresid,map.values)

  #class(vss.results) <- c("psych", "vss")
  if(truemodel){
  return(vss.stats[1,])
    }else{return(vss.stats)}
  #return(f)
}

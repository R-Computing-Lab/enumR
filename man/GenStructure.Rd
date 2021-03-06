% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenStructure.R
\name{GenStructure}
\alias{GenStructure}
\title{Generate Effects matrix}
\usage{
GenStructure(nfactors = 5, rfactors = 0, r_norm = FALSE,
  r_norm_sd = 0.015)
}
\arguments{
\item{nfactors}{Number of factors}

\item{rfactors}{Intercorrelations between factors. Default is 0}

\item{r_norm}{If TRUE, generate normally distributed factor loadings with  mean \code{loading} and sd \code{r_norm_sd}. Does not check if correlation matrix is positive definite.}

\item{r_norm_sd}{If \code{r_norm} is TRUE, standard deviation for normally distributed intercorrelations.}
}
\value{
\code{matrix} of factors and their intercorrelations.
}
\description{
\code{GenStructure} returns an effects matrix, based on intercorrelations and number of factors
}

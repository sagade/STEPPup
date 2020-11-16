\name{CutPointCV}
\alias{CutPointCV}
\title{Helper function for parallel processing to compute prediction error measure for different cuts of the variable}
\usage{
CutPointCV(c, response, var, treat, covars, data, K, repeats, seed)
}
\arguments{
  \item{c}{cutpoint for variable}

  \item{response}{the response, see
  \code{link{OptimalCutPoint}}

  \item{var}{the variable, see \code{link{OptimalCutPoint}}

  \item{treat}{treatment, see \code{link{OptimalCutPoint}}

  \item{covars}{covariates, see
  \code{link{OptimalCutPoint}}

  \item{data}{see \code{link{OptimalCutPoint}}

  \item{K}{see \code{link{OptimalCutPoint}}

  \item{repeats}{see \code{link{OptimalCutPoint}}

  \item{seed}{see \code{link{OptimalCutPoint}}
}
\description{
Helper function for parallel processing to compute
prediction error measure for different cuts of the variable
}

